mod external_data_readers;

pub use external_data_readers::ExternalDataReaders;

use std::{error, fmt, io};

use bytes::Buf;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{self as sam, alignment::record_buf::QualityScores};

use crate::{
    container::block,
    data_container::{
        compression_header::{
            data_series_encoding_map::DataSeries, preservation_map::tag_ids_dictionary,
        },
        CompressionHeader, ReferenceSequenceContext,
    },
    io::BitReader,
    record::{
        feature::{self, substitution},
        Feature, Flags, NextMateFlags,
    },
    Record,
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReadRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_ids_dictionary::Key),
}

impl error::Error for ReadRecordError {}

impl fmt::Display for ReadRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDataSeriesEncoding(data_series) => {
                write!(f, "missing data series encoding: {data_series:?}")
            }
            Self::MissingTagEncoding(key) => write!(f, "missing tag encoding: {key:?}"),
        }
    }
}

pub struct Reader<'a, CDR, EDR>
where
    CDR: Buf,
    EDR: Buf,
{
    compression_header: &'a CompressionHeader,
    core_data_reader: BitReader<CDR>,
    external_data_readers: ExternalDataReaders<EDR>,
    reference_sequence_context: ReferenceSequenceContext,
    prev_alignment_start: Option<Position>,
}

impl<'a, CDR, EDR> Reader<'a, CDR, EDR>
where
    CDR: Buf,
    EDR: Buf,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_reader: BitReader<CDR>,
        external_data_readers: ExternalDataReaders<EDR>,
        reference_sequence_context: ReferenceSequenceContext,
    ) -> Self {
        let initial_alignment_start = match reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.alignment_start()),
            _ => None,
        };

        Self {
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_context,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<()> {
        let bam_bit_flags = self.read_bam_bit_flags()?;
        record.bam_bit_flags = bam_bit_flags;

        let cram_bit_flags = self.read_cram_bit_flags()?;
        record.cram_bit_flags = cram_bit_flags;

        let read_length = self.read_positional_data(record)?;
        self.read_read_names(record)?;
        self.read_mate_data(record, bam_bit_flags, cram_bit_flags)?;

        record.tags = self.read_tag_data()?;

        if bam_bit_flags.is_unmapped() {
            self.read_unmapped_read(record, cram_bit_flags, read_length)?;
        } else {
            self.read_mapped_read(record, cram_bit_flags, read_length)?;
        }

        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    fn read_bam_bit_flags(&mut self) -> io::Result<sam::record::Flags> {
        self.compression_header
            .data_series_encoding_map()
            .bam_bit_flags_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::record::Flags::from)
    }

    fn read_cram_bit_flags(&mut self) -> io::Result<Flags> {
        self.compression_header
            .data_series_encoding_map()
            .cram_bit_flags_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Flags::from)
    }

    fn read_positional_data(&mut self, record: &mut Record) -> io::Result<usize> {
        record.reference_sequence_id = match self.reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.reference_sequence_id()),
            ReferenceSequenceContext::None => None,
            ReferenceSequenceContext::Many => self.read_reference_id()?,
        };

        let read_length = self.read_read_length()?;
        record.read_length = read_length;

        record.alignment_start = self.read_alignment_start()?;
        record.read_group_id = self.read_read_group()?;

        Ok(read_length)
    }

    fn read_reference_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encoding_map()
            .reference_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceId),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                UNMAPPED => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_read_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .read_lengths_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_alignment_start(&mut self) -> io::Result<Option<Position>> {
        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let alignment_start_or_delta = self
            .compression_header
            .data_series_encoding_map()
            .in_seq_positions_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

        let alignment_start = if ap_data_series_delta {
            let prev_alignment_start = i32::try_from(
                self.prev_alignment_start
                    .map(usize::from)
                    .unwrap_or_default(),
            )
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            prev_alignment_start + alignment_start_or_delta
        } else {
            alignment_start_or_delta
        };

        usize::try_from(alignment_start)
            .map(Position::new)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    fn read_read_group(&mut self) -> io::Result<Option<usize>> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i32 = -1;

        self.compression_header
            .data_series_encoding_map()
            .read_groups_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                MISSING => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_read_names(&mut self, record: &mut Record) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();

        // Missing read names are generated when resolving mates.
        if preservation_map.read_names_included() {
            record.name = self.read_read_name()?;
        }

        Ok(())
    }

    fn read_read_name(&mut self) -> io::Result<Option<sam::alignment::record_buf::Name>> {
        const MISSING: &[u8] = &[b'*', 0x00];

        let buf = self
            .compression_header
            .data_series_encoding_map()
            .read_names_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadNames),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

        let name = match &buf[..] {
            MISSING => None,
            _ => Some(sam::alignment::record_buf::Name::from(buf)),
        };

        Ok(name)
    }

    fn read_mate_data(
        &mut self,
        record: &mut Record,
        mut bam_flags: sam::record::Flags,
        flags: Flags,
    ) -> io::Result<()> {
        if flags.is_detached() {
            let next_mate_bit_flags = self.read_next_mate_bit_flags()?;
            record.next_mate_bit_flags = next_mate_bit_flags;

            if next_mate_bit_flags.is_on_negative_strand() {
                bam_flags |= sam::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if next_mate_bit_flags.is_unmapped() {
                bam_flags |= sam::record::Flags::MATE_UNMAPPED;
            }

            record.bam_bit_flags = bam_flags;

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                record.name = self.read_read_name()?;
            }

            record.next_fragment_reference_sequence_id =
                self.read_next_fragment_reference_sequence_id()?;

            record.next_mate_alignment_start = self.read_next_mate_alignment_start()?;
            record.template_size = self.read_template_size()?;
        } else if flags.has_mate_downstream() {
            record.distance_to_next_fragment = self.read_distance_to_next_fragment().map(Some)?;
        }

        Ok(())
    }

    fn read_next_mate_bit_flags(&mut self) -> io::Result<NextMateFlags> {
        self.compression_header
            .data_series_encoding_map()
            .next_mate_bit_flags_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateBitFlags),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(NextMateFlags::from)
    }

    fn read_next_fragment_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encoding_map()
            .next_fragment_reference_sequence_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(
                        DataSeries::NextFragmentReferenceSequenceId,
                    ),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| match id {
                UNMAPPED => Ok(None),
                _ => usize::try_from(id)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_next_mate_alignment_start(&mut self) -> io::Result<Option<Position>> {
        self.compression_header
            .data_series_encoding_map()
            .next_mate_alignment_start_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateAlignmentStart),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Position::new)
    }

    fn read_template_size(&mut self) -> io::Result<i32> {
        self.compression_header
            .data_series_encoding_map()
            .template_size_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::TemplateSize),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_distance_to_next_fragment(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .distance_to_next_fragment_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DistanceToNextFragment),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_tag_data(&mut self) -> io::Result<sam::alignment::record_buf::Data> {
        use bam::record::codec::decoder::data::field::get_value;

        let tag_line = self.read_tag_line()?;

        let tag_keys = self
            .compression_header
            .preservation_map()
            .tag_ids_dictionary()
            .get(tag_line)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid tag line"))?;

        let tag_encoding_map = self.compression_header.tag_encoding_map();

        let mut fields = Vec::with_capacity(tag_keys.len());

        for &key in tag_keys {
            let id = block::ContentId::from(key);

            let encoding = tag_encoding_map.get(&id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingTagEncoding(key),
                )
            })?;

            let data =
                encoding.decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

            let mut data_reader = &data[..];
            let value = get_value(&mut data_reader, key.ty())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let field = (key.tag(), value);
            fields.push(field);
        }

        Ok(fields.into_iter().collect())
    }

    fn read_tag_line(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .tag_ids_encoding()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapped_read(
        &mut self,
        record: &mut Record,
        flags: Flags,
        read_length: usize,
    ) -> io::Result<()> {
        let feature_count = self.read_number_of_read_features()?;

        let mut prev_position = 0;

        for _ in 0..feature_count {
            let feature = self.read_feature(prev_position)?;
            prev_position = usize::from(feature.position());
            record.features.push(feature);
        }

        record.mapping_quality = self.read_mapping_quality()?;

        if flags.are_quality_scores_stored_as_array() {
            record.quality_scores = self.read_quality_scores_stored_as_array(read_length)?;
        }

        Ok(())
    }

    fn read_number_of_read_features(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .number_of_read_features_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NumberOfReadFeatures),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature(&mut self, prev_position: usize) -> io::Result<Feature> {
        use feature::Code;

        let code = self.read_feature_code()?;

        let delta = self.read_feature_position()?;
        let position = Position::try_from(prev_position + delta)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match code {
            Code::Bases => {
                let bases = self.read_stretches_of_bases()?;
                Ok(Feature::Bases(position, bases))
            }
            Code::Scores => {
                let quality_scores = self.read_stretches_of_quality_scores()?;
                Ok(Feature::Scores(position, quality_scores))
            }
            Code::ReadBase => {
                let base = self.read_base()?;
                let quality_score = self.read_quality_score()?;
                Ok(Feature::ReadBase(position, base, quality_score))
            }
            Code::Substitution => {
                let code = self.read_base_substitution_code()?;
                Ok(Feature::Substitution(position, code))
            }
            Code::Insertion => {
                let bases = self.read_insertion()?;
                Ok(Feature::Insertion(position, bases))
            }
            Code::Deletion => {
                let len = self.read_deletion_length()?;
                Ok(Feature::Deletion(position, len))
            }
            Code::InsertBase => {
                let base = self.read_base()?;
                Ok(Feature::InsertBase(position, base))
            }
            Code::QualityScore => {
                let score = self.read_quality_score()?;
                Ok(Feature::QualityScore(position, score))
            }
            Code::ReferenceSkip => {
                let len = self.read_reference_skip_length()?;
                Ok(Feature::ReferenceSkip(position, len))
            }
            Code::SoftClip => {
                let bases = self.read_soft_clip()?;
                Ok(Feature::SoftClip(position, bases))
            }
            Code::Padding => {
                let len = self.read_padding()?;
                Ok(Feature::Padding(position, len))
            }
            Code::HardClip => {
                let len = self.read_hard_clip()?;
                Ok(Feature::HardClip(position, len))
            }
        }
    }

    fn read_feature_code(&mut self) -> io::Result<feature::Code> {
        self.compression_header
            .data_series_encoding_map()
            .read_features_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadFeaturesCodes),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| {
                feature::Code::try_from(id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature_position(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .in_read_positions_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::InReadPositions),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_stretches_of_bases(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encoding_map()
            .stretches_of_bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::StretchesOfBases),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_stretches_of_quality_scores(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encoding_map()
            .stretches_of_quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(
                        DataSeries::StretchesOfQualityScores,
                    ),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encoding_map()
            .bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Bases),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_quality_score(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encoding_map()
            .quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::QualityScores),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base_substitution_code(&mut self) -> io::Result<substitution::Value> {
        self.compression_header
            .data_series_encoding_map()
            .base_substitution_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::BaseSubstitutionCodes),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .map(substitution::Value::Code)
    }

    fn read_insertion(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encoding_map()
            .insertion_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Insertion),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_deletion_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .deletion_lengths_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DeletionLengths),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_reference_skip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .reference_skip_length_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSkipLength),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_soft_clip(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encoding_map()
            .soft_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::SoftClip),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_padding(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .padding_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Padding),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_hard_clip(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encoding_map()
            .hard_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::HardClip),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapping_quality(&mut self) -> io::Result<Option<sam::record::MappingQuality>> {
        self.compression_header
            .data_series_encoding_map()
            .mapping_qualities_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MappingQualities),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::record::MappingQuality::new)
    }

    fn read_unmapped_read(
        &mut self,
        record: &mut Record,
        flags: Flags,
        read_length: usize,
    ) -> io::Result<()> {
        record.bases.as_mut().reserve(read_length);

        for _ in 0..read_length {
            let base = self.read_base()?;
            record.bases.as_mut().push(base);
        }

        if flags.are_quality_scores_stored_as_array() {
            record.quality_scores = self.read_quality_scores_stored_as_array(read_length)?;
        }

        Ok(())
    }

    fn read_quality_scores_stored_as_array(
        &mut self,
        read_length: usize,
    ) -> io::Result<QualityScores> {
        const MISSING: u8 = 0xff;

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::QualityScores),
                )
            })?;

        let mut buf = vec![0; read_length];

        encoding.get().decode_exact(
            &mut self.core_data_reader,
            &mut self.external_data_readers,
            &mut buf,
        )?;

        if buf.iter().all(|&n| n == MISSING) {
            buf.clear();
        }

        Ok(QualityScores::from(buf))
    }
}
