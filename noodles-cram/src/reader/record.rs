use std::{
    collections::HashMap,
    convert::TryFrom,
    error, fmt,
    io::{self, BufRead, Read},
};

use byteorder::ReadBytesExt;
use noodles_bam as bam;
use noodles_sam as sam;

use super::num::read_itf8;
use crate::{
    container::ReferenceSequenceId,
    data_container::{
        compression_header::{data_series_encoding_map::DataSeries, encoding::Encoding},
        CompressionHeader,
    },
    huffman::CanonicalHuffmanDecoder,
    num::Itf8,
    record::{self, feature, tag, Feature, ReadGroupId, Tag},
    BitReader, Record,
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReadRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag::Key),
    MissingExternalBlock(i32),
}

impl error::Error for ReadRecordError {}

impl fmt::Display for ReadRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDataSeriesEncoding(data_series) => {
                write!(f, "missing data series encoding: {:?}", data_series)
            }
            Self::MissingTagEncoding(key) => write!(f, "missing tag encoding: {:?}", key),
            Self::MissingExternalBlock(block_content_id) => {
                write!(f, "missing external block: {}", block_content_id)
            }
        }
    }
}

pub struct Reader<'a, R, S>
where
    R: Read,
    S: BufRead,
{
    compression_header: &'a CompressionHeader,
    core_data_reader: BitReader<R>,
    external_data_readers: HashMap<Itf8, S>,
    reference_sequence_id: ReferenceSequenceId,
    prev_alignment_start: Itf8,
}

impl<'a, R, S> Reader<'a, R, S>
where
    R: Read,
    S: BufRead,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_reader: BitReader<R>,
        external_data_readers: HashMap<Itf8, S>,
        reference_sequence_id: ReferenceSequenceId,
        initial_alignment_start: Itf8,
    ) -> Self {
        Self {
            compression_header,
            core_data_reader,
            external_data_readers,
            reference_sequence_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<()> {
        record.bam_bit_flags = self
            .read_bam_bit_flags()
            .and_then(|n| {
                u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(sam::record::Flags::from)?;

        record.cram_bit_flags = self
            .read_cram_bit_flags()
            .and_then(|n| {
                u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
            .map(record::Flags::from)?;

        self.read_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        // Missing read names are generated when resolving mates.
        if preservation_map.read_names_included() {
            record.read_name = self.read_read_name()?;
        }

        self.read_mate_data(record)?;
        self.read_tag_data(record)?;

        if record.bam_flags().is_unmapped() {
            self.read_unmapped_read(record)?;
        } else {
            self.read_mapped_read(record)?;
        }

        self.prev_alignment_start = record.alignment_start();

        Ok(())
    }

    fn read_bam_bit_flags(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .bam_bit_flags_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_cram_bit_flags(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .cram_bit_flags_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_positional_data(&mut self, record: &mut Record) -> io::Result<()> {
        let reference_id = if self.reference_sequence_id.is_many() {
            self.read_reference_id()?
        } else {
            i32::from(self.reference_sequence_id)
        };

        record.reference_sequence_id =
            if reference_id == bam::record::reference_sequence_id::UNMAPPED {
                None
            } else {
                bam::record::ReferenceSequenceId::try_from(reference_id)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            };

        record.read_length = self.read_read_length()?;

        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let alignment_start = self.read_alignment_start()?;

        record.alignment_start = if ap_data_series_delta {
            self.prev_alignment_start + alignment_start
        } else {
            alignment_start
        };

        record.read_group = self.read_read_group().map(ReadGroupId::from)?;

        Ok(())
    }

    fn read_reference_id(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .reference_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceId),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_read_length(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_lengths_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_alignment_start(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .in_seq_positions_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_read_group(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_groups_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_read_name(&mut self) -> io::Result<Vec<u8>> {
        self.compression_header
            .data_series_encoding_map()
            .read_names_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadNames),
                )
            })
            .and_then(|encoding| {
                decode_byte_array(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                    None,
                )
            })
    }

    fn read_mate_data(&mut self, record: &mut Record) -> io::Result<()> {
        let flags = record.flags();

        if flags.is_detached() {
            record.next_mate_bit_flags = self
                .read_next_mate_bit_flags()
                .and_then(|n| {
                    u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                })
                .map(record::NextMateFlags::from)?;
            let next_mate_flags = record.next_mate_flags();

            if next_mate_flags.is_on_negative_strand() {
                record.bam_bit_flags |= sam::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if next_mate_flags.is_unmapped() {
                record.bam_bit_flags |= sam::record::Flags::MATE_UNMAPPED;
            }

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                record.read_name = self.read_read_name()?;
            }

            record.next_fragment_reference_sequence_id = self
                .read_next_fragment_reference_sequence_id()
                .and_then(|id| {
                    if id == bam::record::reference_sequence_id::UNMAPPED {
                        Ok(None)
                    } else {
                        bam::record::ReferenceSequenceId::try_from(id)
                            .map(Some)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    }
                })?;

            record.next_mate_alignment_start = self.read_next_mate_alignment_start()?;
            record.template_size = self.read_template_size()?;
        } else if flags.has_mate_downstream() {
            record.distance_to_next_fragment = self.read_distance_to_next_fragment()?;
        }

        Ok(())
    }

    fn read_next_mate_bit_flags(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .next_mate_bit_flags_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateBitFlags),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_next_fragment_reference_sequence_id(&mut self) -> io::Result<Itf8> {
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
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_next_mate_alignment_start(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .next_mate_alignment_start_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateAlignmentStart),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_template_size(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .template_size_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::TemplateSize),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_distance_to_next_fragment(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .distance_to_next_fragment_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DistanceToNextFragment),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_tag_data(&mut self, record: &mut Record) -> io::Result<()> {
        let tag_line = self.read_tag_line()?;

        let preservation_map = self.compression_header.preservation_map();
        let tag_ids_dictionary = preservation_map.tag_ids_dictionary();
        let tag_keys = &tag_ids_dictionary[tag_line as usize];

        let tag_encoding_map = self.compression_header.tag_encoding_map();

        record.tags.clear();

        for key in tag_keys {
            let id = key.id();
            let encoding = tag_encoding_map.get(&id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingTagEncoding(*key),
                )
            })?;

            let data = decode_byte_array(
                encoding,
                &mut self.core_data_reader,
                &mut self.external_data_readers,
                None,
            )?;

            let mut data_reader = bam::record::data::Reader::new(&data[..]);
            let value = data_reader.read_value_type(key.ty())?;

            let tag = Tag::new(*key, value);
            record.add_tag(tag);
        }

        Ok(())
    }

    fn read_tag_line(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .tag_ids_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
    }

    fn read_mapped_read(&mut self, record: &mut Record) -> io::Result<()> {
        let read_features_len = self.read_number_of_read_features()?;

        record.features.clear();

        let mut prev_position = 0;

        for _ in 0..read_features_len {
            let feature = self.read_feature(prev_position)?;
            prev_position = feature.position();
            record.add_feature(feature);
        }

        record.mapping_quality = self
            .read_mapping_quality()
            .map(|n| sam::record::MappingQuality::from(n as u8))?;

        let flags = record.flags();

        if flags.are_quality_scores_stored_as_array() {
            let read_len = record.read_length();
            let mut scores = Vec::with_capacity(read_len as usize);

            for _ in 0..read_len {
                let score = self.read_quality_score()?;
                scores.push(score);
            }

            record.quality_scores = scores;
        }

        Ok(())
    }

    fn read_number_of_read_features(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .number_of_read_features_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NumberOfReadFeatures),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_feature(&mut self, prev_position: i32) -> io::Result<Feature> {
        let code = self
            .read_feature_code()
            .map(|id| id as u8 as char)
            .and_then(|id| {
                feature::Code::try_from(id)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

        let position = prev_position + self.read_feature_position()?;

        match code {
            feature::Code::Bases => {
                let bases = self.read_stretches_of_bases()?;
                Ok(Feature::Bases(position, bases))
            }
            feature::Code::Scores => {
                let quality_scores = self.read_stretches_of_quality_scores()?;
                Ok(Feature::Scores(position, quality_scores))
            }
            feature::Code::ReadBase => {
                let base = self.read_base()?;
                let quality_score = self.read_quality_score()?;
                Ok(Feature::ReadBase(position, base, quality_score))
            }
            feature::Code::Substitution => {
                let code = self.read_base_substitution_code()?;
                Ok(Feature::Substitution(position, code))
            }
            feature::Code::Insertion => {
                let bases = self.read_insertion()?;
                Ok(Feature::Insertion(position, bases))
            }
            feature::Code::Deletion => {
                let len = self.read_deletion_length()?;
                Ok(Feature::Deletion(position, len))
            }
            feature::Code::InsertBase => {
                let base = self.read_base()?;
                Ok(Feature::InsertBase(position, base))
            }
            feature::Code::QualityScore => {
                let score = self.read_quality_score()?;
                Ok(Feature::QualityScore(position, score))
            }
            feature::Code::ReferenceSkip => {
                let len = self.read_reference_skip_length()?;
                Ok(Feature::ReferenceSkip(position, len))
            }
            feature::Code::SoftClip => {
                let bases = self.read_soft_clip()?;
                Ok(Feature::SoftClip(position, bases))
            }
            feature::Code::Padding => {
                let len = self.read_padding()?;
                Ok(Feature::Padding(position, len))
            }
            feature::Code::HardClip => {
                let len = self.read_hard_clip()?;
                Ok(Feature::HardClip(position, len))
            }
        }
    }

    fn read_feature_code(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .read_features_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadFeaturesCodes),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_feature_position(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .in_read_positions_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::InReadPositions),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
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
            })
            .and_then(|encoding| {
                decode_byte_array(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                    None,
                )
            })
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
            })
            .and_then(|encoding| {
                decode_byte_array(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                    None,
                )
            })
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
            })
            .and_then(|encoding| {
                decode_byte(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
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
            })
            .and_then(|encoding| {
                decode_byte(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_base_substitution_code(&mut self) -> io::Result<u8> {
        self.compression_header
            .data_series_encoding_map()
            .base_substitution_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::BaseSubstitutionCodes),
                )
            })
            .and_then(|encoding| {
                decode_byte(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
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
            })
            .and_then(|encoding| {
                decode_byte_array(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                    None,
                )
            })
    }

    fn read_deletion_length(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .deletion_lengths_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DeletionLengths),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_reference_skip_length(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .reference_skip_length_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSkipLength),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
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
            })
            .and_then(|encoding| {
                decode_byte_array(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                    None,
                )
            })
    }

    fn read_padding(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .padding_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Padding),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_hard_clip(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .hard_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::HardClip),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_mapping_quality(&mut self) -> io::Result<Itf8> {
        self.compression_header
            .data_series_encoding_map()
            .mapping_qualities_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MappingQualities),
                )
            })
            .and_then(|encoding| {
                decode_itf8(
                    encoding,
                    &mut self.core_data_reader,
                    &mut self.external_data_readers,
                )
            })
    }

    fn read_unmapped_read(&mut self, record: &mut Record) -> io::Result<()> {
        record.bases.clear();

        for _ in 0..record.read_length() {
            let base = self.read_base()?;
            record.bases.push(base);
        }

        let flags = record.flags();

        if flags.are_quality_scores_stored_as_array() {
            let read_len = record.read_length();
            let mut scores = Vec::with_capacity(read_len as usize);

            for _ in 0..read_len {
                let score = self.read_quality_score()?;
                scores.push(score);
            }

            record.quality_scores = scores;
        }

        Ok(())
    }
}

fn decode_byte<R, S>(
    encoding: &Encoding,
    core_data_reader: &mut BitReader<R>,
    external_data_readers: &mut HashMap<Itf8, S>,
) -> io::Result<u8>
where
    R: Read,
    S: Read,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let reader = external_data_readers
                .get_mut(block_content_id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingExternalBlock(*block_content_id),
                    )
                })?;

            reader.read_u8()
        }
        Encoding::Huffman(alphabet, bit_lens) => {
            let decoder = CanonicalHuffmanDecoder::new(alphabet, bit_lens);
            decoder.read(core_data_reader).map(|i| i as u8)
        }
        Encoding::Beta(offset, len) => core_data_reader
            .read_u32(*len as usize)
            .map(|i| (i as i32 - offset) as u8),
        _ => todo!("decode_byte: {:?}", encoding),
    }
}

fn decode_itf8<R, S>(
    encoding: &Encoding,
    core_data_reader: &mut BitReader<R>,
    external_data_readers: &mut HashMap<Itf8, S>,
) -> io::Result<Itf8>
where
    R: Read,
    S: Read,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let reader = external_data_readers
                .get_mut(block_content_id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingExternalBlock(*block_content_id),
                    )
                })?;

            read_itf8(reader)
        }
        Encoding::Huffman(alphabet, bit_lens) => {
            let decoder = CanonicalHuffmanDecoder::new(alphabet, bit_lens);
            decoder.read(core_data_reader)
        }
        Encoding::Beta(offset, len) => core_data_reader
            .read_u32(*len as usize)
            .map(|i| (i as i32 - offset)),
        _ => todo!("decode_itf8: {:?}", encoding),
    }
}

fn decode_byte_array<R, S>(
    encoding: &Encoding,
    core_data_reader: &mut BitReader<R>,
    external_data_readers: &mut HashMap<Itf8, S>,
    buf: Option<Vec<u8>>,
) -> io::Result<Vec<u8>>
where
    R: Read,
    S: BufRead,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let reader = external_data_readers
                .get_mut(block_content_id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingExternalBlock(*block_content_id),
                    )
                })?;

            let mut buf = buf.expect("missing buf");
            reader.read_exact(&mut buf)?;

            Ok(buf)
        }
        Encoding::ByteArrayLen(len_encoding, value_encoding) => {
            let len = decode_itf8(len_encoding, core_data_reader, external_data_readers)?;

            let buf = vec![0; len as usize];
            let value = decode_byte_array(
                value_encoding,
                core_data_reader,
                external_data_readers,
                Some(buf),
            )?;

            Ok(value)
        }
        Encoding::ByteArrayStop(stop_byte, block_content_id) => {
            let reader = external_data_readers
                .get_mut(block_content_id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        ReadRecordError::MissingExternalBlock(*block_content_id),
                    )
                })?;

            let mut buf = Vec::new();
            reader.read_until(*stop_byte, &mut buf)?;

            // Remove stop byte.
            buf.pop();

            Ok(buf)
        }
        _ => todo!("decode_byte_array: {:?}", encoding),
    }
}
