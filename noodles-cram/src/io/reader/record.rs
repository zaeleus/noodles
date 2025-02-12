mod external_data_readers;

pub use external_data_readers::ExternalDataReaders;

use std::{error, fmt, io};

use bstr::BString;
use bytes::Buf;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{self as sam, alignment::record_buf::QualityScores};

use crate::{
    container::{
        block,
        compression_header::{data_series_encodings::DataSeries, preservation_map::tag_sets},
        CompressionHeader, ReferenceSequenceContext,
    },
    io::BitReader,
    record::{feature, Feature, Flags, MateFlags},
    Record,
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReadRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_sets::Key),
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
        record.bam_flags = self.read_bam_flags()?;
        record.cram_flags = self.read_cram_flags()?;

        self.read_positions(record)?;
        self.read_names(record)?;
        self.read_mate(record)?;

        record.tags = self.read_tags()?;

        if record.bam_flags().is_unmapped() {
            self.read_unmapped_read(record)?;
        } else {
            self.read_mapped_read(record)?;
        }

        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    fn read_bam_flags(&mut self) -> io::Result<sam::alignment::record::Flags> {
        self.compression_header
            .data_series_encodings()
            .bam_flags()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::alignment::record::Flags::from)
    }

    fn read_cram_flags(&mut self) -> io::Result<Flags> {
        self.compression_header
            .data_series_encodings()
            .cram_flags()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Flags::from)
    }

    fn read_positions(&mut self, record: &mut Record) -> io::Result<()> {
        record.reference_sequence_id = match self.reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.reference_sequence_id()),
            ReferenceSequenceContext::None => None,
            ReferenceSequenceContext::Many => self.read_reference_sequence_id()?,
        };

        let read_length = self.read_read_length()?;
        record.read_length = read_length;

        record.alignment_start = self.read_alignment_start()?;
        record.read_group_id = self.read_read_group_id()?;

        Ok(())
    }

    fn read_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .reference_sequence_ids()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSequenceIds),
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
            .data_series_encodings()
            .read_lengths()
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
            .data_series_encodings()
            .alignment_starts()
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

    fn read_read_group_id(&mut self) -> io::Result<Option<usize>> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .read_group_ids()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| match n {
                MISSING => Ok(None),
                _ => usize::try_from(n)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
            })
    }

    fn read_names(&mut self, record: &mut Record) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();

        // Missing read names are generated when resolving mates.
        if preservation_map.read_names_included() {
            record.name = self.read_name()?;
        }

        Ok(())
    }

    fn read_name(&mut self) -> io::Result<Option<BString>> {
        const MISSING: &[u8] = &[b'*', 0x00];

        let buf = self
            .compression_header
            .data_series_encodings()
            .names()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Names),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

        let name = match &buf[..] {
            MISSING => None,
            _ => Some(BString::from(buf)),
        };

        Ok(name)
    }

    fn read_mate(&mut self, record: &mut Record) -> io::Result<()> {
        if record.cram_flags().is_detached() {
            record.mate_flags = self.read_mate_flags()?;

            if record.next_mate_flags().is_on_negative_strand() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if record.next_mate_flags().is_unmapped() {
                record.bam_flags |= sam::alignment::record::Flags::MATE_UNMAPPED;
            }

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                record.name = self.read_name()?;
            }

            record.mate_reference_sequence_id = self.read_mate_reference_sequence_id()?;

            record.mate_alignment_start = self.read_mate_alignment_start()?;
            record.template_length = self.read_template_length()?;
        } else if record.cram_flags().has_mate_downstream() {
            record.distance_to_mate = self.read_mate_distance().map(Some)?;
        }

        Ok(())
    }

    fn read_mate_flags(&mut self) -> io::Result<MateFlags> {
        self.compression_header
            .data_series_encodings()
            .mate_flags()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MateFlags),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(MateFlags::from)
    }

    fn read_mate_reference_sequence_id(&mut self) -> io::Result<Option<usize>> {
        const UNMAPPED: i32 = -1;

        self.compression_header
            .data_series_encodings()
            .mate_reference_sequence_ids()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MateReferenceSequenceId),
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

    fn read_mate_alignment_start(&mut self) -> io::Result<Option<Position>> {
        self.compression_header
            .data_series_encodings()
            .mate_alignment_starts()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MateAlignmentStart),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(Position::new)
    }

    fn read_template_length(&mut self) -> io::Result<i32> {
        self.compression_header
            .data_series_encodings()
            .template_lengths()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::TemplateLengths),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_mate_distance(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .mate_distances()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MateDistances),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_tags(&mut self) -> io::Result<sam::alignment::record_buf::Data> {
        use bam::record::codec::decoder::data::field::read_value;

        let tag_set_id = self.read_tag_set_id()?;

        let tag_keys = self
            .compression_header
            .preservation_map()
            .tag_sets()
            .get(tag_set_id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid tag line"))?;

        let tag_encodings = self.compression_header.tag_encodings();

        let mut fields = Vec::with_capacity(tag_keys.len());

        for &key in tag_keys {
            let id = block::ContentId::from(key);

            let encoding = tag_encodings.get(&id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingTagEncoding(key),
                )
            })?;

            let data =
                encoding.decode(&mut self.core_data_reader, &mut self.external_data_readers)?;

            let mut data_reader = &data[..];
            let value = read_value(&mut data_reader, key.ty())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let field = (key.tag(), value);
            fields.push(field);
        }

        Ok(fields.into_iter().collect())
    }

    fn read_tag_set_id(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .tag_set_ids()
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapped_read(&mut self, record: &mut Record) -> io::Result<()> {
        let feature_count = self.read_feature_count()?;

        let mut prev_position = 0;

        for _ in 0..feature_count {
            let feature = self.read_feature(prev_position)?;
            prev_position = usize::from(feature.position());
            record.features.push(feature);
        }

        record.mapping_quality = self.read_mapping_quality()?;

        if record.cram_flags().are_quality_scores_stored_as_array() {
            record.quality_scores =
                self.read_quality_scores_stored_as_array(record.read_length())?;
        }

        Ok(())
    }

    fn read_feature_count(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_counts()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::FeatureCounts),
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

        let delta = self.read_feature_position_delta()?;
        let position = Position::try_from(prev_position + delta)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        match code {
            Code::Bases => {
                let bases = self.read_stretches_of_bases()?;
                Ok(Feature::Bases { position, bases })
            }
            Code::Scores => {
                let quality_scores = self.read_stretches_of_quality_scores()?;

                Ok(Feature::Scores {
                    position,
                    quality_scores,
                })
            }
            Code::ReadBase => {
                let base = self.read_base()?;
                let quality_score = self.read_quality_score()?;

                Ok(Feature::ReadBase {
                    position,
                    base,
                    quality_score,
                })
            }
            Code::Substitution => {
                let code = self.read_base_substitution_code()?;
                Ok(Feature::Substitution { position, code })
            }
            Code::Insertion => {
                let bases = self.read_insertion_bases()?;
                Ok(Feature::Insertion { position, bases })
            }
            Code::Deletion => {
                let len = self.read_deletion_length()?;
                Ok(Feature::Deletion { position, len })
            }
            Code::InsertBase => {
                let base = self.read_base()?;
                Ok(Feature::InsertBase { position, base })
            }
            Code::QualityScore => {
                let quality_score = self.read_quality_score()?;

                Ok(Feature::QualityScore {
                    position,
                    quality_score,
                })
            }
            Code::ReferenceSkip => {
                let len = self.read_reference_skip_length()?;
                Ok(Feature::ReferenceSkip { position, len })
            }
            Code::SoftClip => {
                let bases = self.read_soft_clip_bases()?;
                Ok(Feature::SoftClip { position, bases })
            }
            Code::Padding => {
                let len = self.read_padding_length()?;
                Ok(Feature::Padding { position, len })
            }
            Code::HardClip => {
                let len = self.read_hard_clip_length()?;
                Ok(Feature::HardClip { position, len })
            }
        }
    }

    fn read_feature_code(&mut self) -> io::Result<feature::Code> {
        self.compression_header
            .data_series_encodings()
            .feature_codes()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::FeatureCodes),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|id| {
                feature::Code::try_from(id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_feature_position_delta(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .feature_position_deltas()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::FeaturePositionDeltas),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_stretches_of_bases(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_bases()
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
            .data_series_encodings()
            .stretches_of_quality_scores()
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
            .data_series_encodings()
            .bases()
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
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::QualityScores),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_base_substitution_code(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encodings()
            .base_substitution_codes()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::BaseSubstitutionCodes),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_insertion_bases(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encodings()
            .insertion_bases()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::InsertionBases),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_deletion_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .deletion_lengths()
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
            .data_series_encodings()
            .reference_skip_lengths()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSkipLengths),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_soft_clip_bases(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encodings()
            .soft_clip_bases()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::SoftClipBases),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
    }

    fn read_padding_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .padding_lengths()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::PaddingLengths),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_hard_clip_length(&mut self) -> io::Result<usize> {
        self.compression_header
            .data_series_encodings()
            .hard_clip_lengths()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::HardClipLengths),
                )
            })?
            .decode(&mut self.core_data_reader, &mut self.external_data_readers)
            .and_then(|n| {
                usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    fn read_mapping_quality(
        &mut self,
    ) -> io::Result<Option<sam::alignment::record::MappingQuality>> {
        self.compression_header
            .data_series_encodings()
            .mapping_qualities()
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
            .map(sam::alignment::record::MappingQuality::new)
    }

    fn read_unmapped_read(&mut self, record: &mut Record) -> io::Result<()> {
        let read_length = record.read_length();

        record.sequence.as_mut().reserve(read_length);

        for _ in 0..read_length {
            let base = self.read_base()?;
            record.sequence.as_mut().push(base);
        }

        if record.cram_flags().are_quality_scores_stored_as_array() {
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
            .data_series_encodings()
            .quality_scores()
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
