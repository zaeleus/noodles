use std::{collections::HashMap, error, fmt, io};

use bstr::BStr;
use noodles_core::Position;
use noodles_sam as sam;

use crate::{
    container::{
        CompressionHeader, ReferenceSequenceContext, block,
        compression_header::{
            data_series_encodings::DataSeries,
            preservation_map::{substitution_matrix::Base, tag_sets},
        },
    },
    io::{
        BitWriter,
        writer::{Record, record::Feature},
    },
    record::{Flags, MateFlags},
};

pub type ExternalDataWriters = HashMap<block::ContentId, Vec<u8>>;

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum WriteRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_sets::Key),
}

impl error::Error for WriteRecordError {}

impl fmt::Display for WriteRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDataSeriesEncoding(data_series) => {
                write!(f, "missing data series encoding: {data_series:?}")
            }
            Self::MissingTagEncoding(key) => write!(f, "missing tag encoding: {key:?}"),
        }
    }
}

pub struct Writer<'a> {
    compression_header: &'a CompressionHeader,
    core_data_writer: &'a mut BitWriter,
    external_data_writers: &'a mut ExternalDataWriters,
    reference_sequence_context: ReferenceSequenceContext,
    prev_alignment_start: Option<Position>,
}

impl<'a> Writer<'a> {
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_writer: &'a mut BitWriter,
        external_data_writers: &'a mut ExternalDataWriters,
        reference_sequence_context: ReferenceSequenceContext,
    ) -> Self {
        let initial_alignment_start = match reference_sequence_context {
            ReferenceSequenceContext::Some(context) => Some(context.alignment_start()),
            _ => None,
        };

        Self {
            compression_header,
            core_data_writer,
            external_data_writers,
            reference_sequence_context,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write_bam_flags(record.bam_flags)?;
        self.write_cram_flags(record.cram_flags)?;

        self.write_positions(record)?;
        self.write_names(record)?;
        self.write_mate(record)?;

        self.write_data(record)?;

        if record.bam_flags.is_unmapped() {
            self.write_unmapped_read(record)?;
        } else {
            self.write_mapped_read(record)?;
        }

        self.prev_alignment_start = record.alignment_start;

        Ok(())
    }

    fn write_bam_flags(&mut self, bam_flags: sam::alignment::record::Flags) -> io::Result<()> {
        let n = i64::from(u16::from(bam_flags));

        self.compression_header
            .data_series_encodings()
            .bam_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BamFlags))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_cram_flags(&mut self, flags: Flags) -> io::Result<()> {
        let n = i64::from(u8::from(flags));

        self.compression_header
            .data_series_encodings()
            .cram_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::CramFlags))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_positions(&mut self, record: &Record) -> io::Result<()> {
        if self.reference_sequence_context.is_many() {
            self.write_reference_sequence_id(record.reference_sequence_id)?;
        }

        self.write_read_length(record.read_length)?;
        self.write_alignment_start(record.alignment_start)?;
        self.write_read_group_id(record.read_group_id)?;

        Ok(())
    }

    fn write_reference_sequence_id(
        &mut self,
        reference_sequence_id: Option<usize>,
    ) -> io::Result<()> {
        const UNMAPPED: i64 = -1;

        let n = if let Some(id) = reference_sequence_id {
            i64::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            UNMAPPED
        };

        self.compression_header
            .data_series_encodings()
            .reference_sequence_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSequenceIds))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_read_length(&mut self, read_length: usize) -> io::Result<()> {
        let n = i64::try_from(read_length)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        self.compression_header
            .data_series_encodings()
            .read_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadLengths))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_alignment_start(&mut self, alignment_start: Option<Position>) -> io::Result<()> {
        const MISSING: i64 = 0;

        fn position_to_i64(position: Position) -> io::Result<i64> {
            i64::try_from(usize::from(position))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        }

        let alignment_starts_are_deltas = self
            .compression_header
            .preservation_map()
            .alignment_starts_are_deltas();

        let alignment_start_or_delta = if alignment_starts_are_deltas {
            let start = alignment_start
                .map(position_to_i64)
                .transpose()?
                .unwrap_or(MISSING);

            let prev_start = self
                .prev_alignment_start
                .map(position_to_i64)
                .transpose()?
                .unwrap_or(MISSING);

            start - prev_start
        } else if let Some(position) = alignment_start {
            position_to_i64(position)?
        } else {
            MISSING
        };

        self.compression_header
            .data_series_encodings()
            .alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::AlignmentStarts))?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                alignment_start_or_delta,
            )
    }

    fn write_read_group_id(&mut self, read_group_id: Option<usize>) -> io::Result<()> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i64 = -1;

        let n = if let Some(id) = read_group_id {
            i64::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            MISSING
        };

        self.compression_header
            .data_series_encodings()
            .read_group_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReadGroupIds))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_names(&mut self, record: &Record) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.records_have_names() {
            let name = record.name.as_ref().map(|s| s.as_ref());
            self.write_name(name)?;
        }

        Ok(())
    }

    fn write_name(&mut self, name: Option<&BStr>) -> io::Result<()> {
        const MISSING: &[u8] = &[b'*', 0x00];

        let buf = name.map(|s| s.as_ref()).unwrap_or(MISSING);

        self.compression_header
            .data_series_encodings()
            .names()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Names))?
            .encode(self.core_data_writer, self.external_data_writers, buf)
    }

    fn write_mate(&mut self, record: &Record) -> io::Result<()> {
        if record.cram_flags.is_detached() {
            self.write_mate_flags(record.mate_flags)?;

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.records_have_names() {
                let name = record.name.as_ref().map(|s| s.as_ref());
                self.write_name(name)?;
            }

            self.write_mate_reference_sequence_id(record.mate_reference_sequence_id)?;
            self.write_mate_alignment_start(record.mate_alignment_start)?;
            self.write_template_length(record.template_length)?;
        } else if let Some(mate_distance) = record.mate_distance {
            self.write_mate_distance(mate_distance)?;
        }

        Ok(())
    }

    fn write_mate_flags(&mut self, mate_flags: MateFlags) -> io::Result<()> {
        let n = i64::from(u8::from(mate_flags));

        self.compression_header
            .data_series_encodings()
            .mate_flags()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateFlags))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mate_reference_sequence_id(
        &mut self,
        mate_reference_sequence_id: Option<usize>,
    ) -> io::Result<()> {
        const UNMAPPED: i64 = -1;

        let n = if let Some(id) = mate_reference_sequence_id {
            i64::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            UNMAPPED
        };

        self.compression_header
            .data_series_encodings()
            .mate_reference_sequence_ids()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::MateReferenceSequenceIds)
            })?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mate_alignment_start(
        &mut self,
        mate_alignment_start: Option<Position>,
    ) -> io::Result<()> {
        const MISSING: i64 = 0;

        let n = if let Some(position) = mate_alignment_start {
            i64::try_from(usize::from(position))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            MISSING
        };

        self.compression_header
            .data_series_encodings()
            .mate_alignment_starts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateAlignmentStarts))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_template_length(&mut self, template_length: i64) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .template_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TemplateLengths))?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                template_length,
            )
    }

    fn write_mate_distance(&mut self, mate_distance: usize) -> io::Result<()> {
        let n = i64::try_from(mate_distance)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .mate_distances()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MateDistances))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_data(&mut self, record: &Record) -> io::Result<()> {
        use noodles_bam::record::codec::encoder::data::field::write_value;

        let tag_set: Vec<_> = record
            .data
            .iter()
            .map(|(tag, value)| tag_sets::Key::new(*tag, value.ty()))
            .collect();

        let tag_set_id = self
            .compression_header
            .preservation_map()
            .tag_sets()
            .iter()
            .position(|set| **set == tag_set)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing tag set"))?;

        self.write_tag_set_id(tag_set_id)?;

        let tag_encodings = self.compression_header.tag_encodings();
        let mut buf = Vec::new();

        for (key, (_, value)) in tag_set.into_iter().zip(&record.data) {
            let block_content_id = block::ContentId::from(key);

            write_value(&mut buf, &value.into())?;

            tag_encodings
                .get(&block_content_id)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        WriteRecordError::MissingTagEncoding(key),
                    )
                })?
                .encode(self.core_data_writer, self.external_data_writers, &buf)?;

            buf.clear();
        }

        Ok(())
    }

    fn write_tag_set_id(&mut self, tag_set_id: usize) -> io::Result<()> {
        let n = i64::try_from(tag_set_id)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .tag_set_ids()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::TagSetIds))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mapped_read(&mut self, record: &Record) -> io::Result<()> {
        self.write_feature_count(record.features.len())?;

        let mut prev_position = 0;

        for feature in &record.features {
            let position = usize::from(feature.position()) - prev_position;
            self.write_feature(feature, position)?;
            prev_position = usize::from(feature.position());
        }

        self.write_mapping_quality(record.mapping_quality)?;

        if record.cram_flags.quality_scores_are_stored_as_array() {
            self.write_quality_scores(&record.quality_scores)?;
        }

        Ok(())
    }

    fn write_feature_count(&mut self, feature_count: usize) -> io::Result<()> {
        let n = i64::try_from(feature_count)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .feature_counts()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCounts))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_feature(&mut self, feature: &Feature, position: usize) -> io::Result<()> {
        self.write_feature_code(feature.code())?;
        self.write_feature_position_delta(position)?;

        match feature {
            Feature::Bases { bases, .. } => self.write_stretches_of_bases(bases)?,
            Feature::Scores { quality_scores, .. } => {
                self.write_stretches_of_quality_scores(quality_scores)?
            }
            Feature::ReadBase {
                base,
                quality_score,
                ..
            } => {
                self.write_base(*base)?;
                self.write_quality_score(*quality_score)?;
            }
            Feature::Substitution {
                reference_base,
                read_base,
                ..
            } => self.write_base_substitution_code(*reference_base, *read_base)?,
            Feature::Insertion { bases, .. } => self.write_insertion_bases(bases)?,
            Feature::Deletion { len, .. } => self.write_deletion_length(*len)?,
            Feature::InsertBase { base, .. } => self.write_base(*base)?,
            Feature::QualityScore { quality_score, .. } => {
                self.write_quality_score(*quality_score)?
            }
            Feature::ReferenceSkip { len, .. } => self.write_reference_skip_length(*len)?,
            Feature::SoftClip { bases, .. } => self.write_soft_clip(bases)?,
            Feature::Padding { len, .. } => self.write_padding_length(*len)?,
            Feature::HardClip { len, .. } => self.write_hard_clip_length(*len)?,
        }

        Ok(())
    }

    fn write_feature_code(&mut self, code: u8) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .feature_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeatureCodes))?
            .encode(self.core_data_writer, self.external_data_writers, code)
    }

    fn write_feature_position_delta(&mut self, delta: usize) -> io::Result<()> {
        let n = i64::try_from(delta).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .feature_position_deltas()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::FeaturePositionDeltas))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_stretches_of_bases(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::StretchesOfBases))?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_stretches_of_quality_scores(&mut self, quality_scores: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .stretches_of_quality_scores()
            .ok_or_else(|| {
                missing_data_series_encoding_error(DataSeries::StretchesOfQualityScores)
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                quality_scores,
            )
    }

    fn write_base(&mut self, base: u8) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?
            .encode(self.core_data_writer, self.external_data_writers, base)
    }

    fn write_quality_score(&mut self, quality_score: u8) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                quality_score,
            )
    }

    fn write_base_substitution_code(
        &mut self,
        reference_base: Base,
        read_base: Base,
    ) -> io::Result<()> {
        let code = self
            .compression_header
            .preservation_map()
            .substitution_matrix()
            .find(reference_base, read_base);

        self.compression_header
            .data_series_encodings()
            .base_substitution_codes()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::BaseSubstitutionCodes))?
            .encode(self.core_data_writer, self.external_data_writers, code)
    }

    fn write_insertion_bases(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .insertion_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::InsertionBases))?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_deletion_length(&mut self, len: usize) -> io::Result<()> {
        let n = i64::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .deletion_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::DeletionLengths))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_reference_skip_length(&mut self, len: usize) -> io::Result<()> {
        let n = i64::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .reference_skip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::ReferenceSkipLengths))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_soft_clip(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .soft_clip_bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::SoftClipBases))?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_padding_length(&mut self, len: usize) -> io::Result<()> {
        let n = i64::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .padding_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::PaddingLengths))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_hard_clip_length(&mut self, len: usize) -> io::Result<()> {
        let n = i64::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        self.compression_header
            .data_series_encodings()
            .hard_clip_lengths()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::HardClipLengths))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mapping_quality(
        &mut self,
        mapping_quality: Option<sam::alignment::record::MappingQuality>,
    ) -> io::Result<()> {
        const MISSING: u8 = 0xff;

        let n = i64::from(mapping_quality.map(u8::from).unwrap_or(MISSING));

        self.compression_header
            .data_series_encodings()
            .mapping_qualities()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::MappingQualities))?
            .encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_unmapped_read(&mut self, record: &Record) -> io::Result<()> {
        self.write_sequence(&record.sequence)?;

        if record.cram_flags.quality_scores_are_stored_as_array() {
            self.write_quality_scores(&record.quality_scores)?;
        }

        Ok(())
    }

    fn write_sequence(&mut self, sequence: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .bases()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::Bases))?
            .get()
            .encode_extend(self.core_data_writer, self.external_data_writers, sequence)
    }

    fn write_quality_scores(&mut self, quality_scores: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encodings()
            .quality_scores()
            .ok_or_else(|| missing_data_series_encoding_error(DataSeries::QualityScores))?
            .get()
            .encode_extend(
                self.core_data_writer,
                self.external_data_writers,
                quality_scores,
            )
    }
}

fn missing_data_series_encoding_error(data_series: DataSeries) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        WriteRecordError::MissingDataSeriesEncoding(data_series),
    )
}
