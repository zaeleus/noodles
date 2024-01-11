use std::{
    collections::HashMap,
    error, fmt,
    io::{self, Write},
};

use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;

use crate::{
    container::block,
    data_container::{
        compression_header::{
            data_series_encoding_map::DataSeries, preservation_map::tag_ids_dictionary,
        },
        CompressionHeader, ReferenceSequenceContext,
    },
    io::BitWriter,
    record::{
        feature::{self, substitution},
        Feature, Flags, NextMateFlags,
    },
    Record,
};

#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum WriteRecordError {
    MissingDataSeriesEncoding(DataSeries),
    MissingTagEncoding(tag_ids_dictionary::Key),
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

pub struct Writer<'a, W, X> {
    compression_header: &'a CompressionHeader,
    core_data_writer: &'a mut BitWriter<W>,
    external_data_writers: &'a mut HashMap<block::ContentId, X>,
    reference_sequence_context: ReferenceSequenceContext,
    prev_alignment_start: Option<Position>,
}

impl<'a, W, X> Writer<'a, W, X>
where
    W: Write,
    X: Write,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_writer: &'a mut BitWriter<W>,
        external_data_writers: &'a mut HashMap<block::ContentId, X>,
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
        self.write_bam_bit_flags(record.bam_flags())?;
        self.write_cram_bit_flags(record.cram_flags())?;

        self.write_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            self.write_read_name(record.name())?;
        }

        self.write_mate_data(record)?;
        self.write_tag_data(record)?;

        if record.bam_flags().is_unmapped() {
            self.write_unmapped_read(record)?;
        } else {
            self.write_mapped_read(record)?;
        }

        self.prev_alignment_start = record.alignment_start();

        Ok(())
    }

    fn write_bam_bit_flags(
        &mut self,
        bam_flags: sam::alignment::record_buf::Flags,
    ) -> io::Result<()> {
        let bam_bit_flags = i32::from(u16::from(bam_flags));

        self.compression_header
            .data_series_encoding_map()
            .bam_bit_flags_encoding()
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                bam_bit_flags,
            )
    }

    fn write_cram_bit_flags(&mut self, flags: Flags) -> io::Result<()> {
        let cram_bit_flags = i32::from(u8::from(flags));

        self.compression_header
            .data_series_encoding_map()
            .cram_bit_flags_encoding()
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                cram_bit_flags,
            )
    }

    fn write_positional_data(&mut self, record: &Record) -> io::Result<()> {
        if self.reference_sequence_context.is_many() {
            self.write_reference_id(record.reference_sequence_id())?;
        }

        self.write_read_length(record.read_length())?;
        self.write_alignment_start(record.alignment_start)?;
        self.write_read_group(record.read_group_id())?;

        Ok(())
    }

    fn write_reference_id(&mut self, reference_sequence_id: Option<usize>) -> io::Result<()> {
        const UNMAPPED: i32 = -1;

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .reference_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceId),
                )
            })?;

        let reference_id = if let Some(id) = reference_sequence_id {
            i32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            UNMAPPED
        };

        encoding.encode(
            self.core_data_writer,
            self.external_data_writers,
            reference_id,
        )
    }

    fn write_read_length(&mut self, read_length: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_lengths_encoding();

        let len = i32::try_from(read_length)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, len)
    }

    fn write_alignment_start(&mut self, alignment_start: Option<Position>) -> io::Result<()> {
        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .in_seq_positions_encoding();

        let alignment_start_or_delta = if ap_data_series_delta {
            match (alignment_start, self.prev_alignment_start) {
                (None, None) => 0,
                (Some(alignment_start), Some(prev_alignment_start)) => {
                    let alignment_start = i32::try_from(usize::from(alignment_start))
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

                    let prev_alignment_start = i32::try_from(usize::from(prev_alignment_start))
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

                    alignment_start - prev_alignment_start
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "invalid alignment start ({:?}) or previous alignment start ({:?})",
                            alignment_start, self.prev_alignment_start
                        ),
                    ));
                }
            }
        } else {
            i32::try_from(alignment_start.map(usize::from).unwrap_or_default())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        };

        encoding.encode(
            self.core_data_writer,
            self.external_data_writers,
            alignment_start_or_delta,
        )
    }

    fn write_read_group(&mut self, read_group_id: Option<usize>) -> io::Result<()> {
        // § 10.2 "CRAM positional data" (2021-10-15): "-1 for no group".
        const MISSING: i32 = -1;

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_groups_encoding();

        let read_group = if let Some(id) = read_group_id {
            i32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
        } else {
            MISSING
        };

        encoding.encode(
            self.core_data_writer,
            self.external_data_writers,
            read_group,
        )
    }

    fn write_read_name(
        &mut self,
        name: Option<&sam::alignment::record_buf::Name>,
    ) -> io::Result<()> {
        const MISSING: &[u8] = &[b'*', 0x00];

        let buf = name.map(|name| name.as_ref()).unwrap_or(MISSING);

        self.compression_header
            .data_series_encoding_map()
            .read_names_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::ReadNames),
                )
            })?
            .encode(self.core_data_writer, self.external_data_writers, buf)
    }

    fn write_mate_data(&mut self, record: &Record) -> io::Result<()> {
        if record.cram_flags().is_detached() {
            self.write_next_mate_bit_flags(record.next_mate_flags())?;

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                self.write_read_name(record.name())?;
            }

            self.write_next_fragment_reference_sequence_id(
                record.next_fragment_reference_sequence_id(),
            )?;

            self.write_next_mate_alignment_start(record.next_mate_alignment_start())?;
            self.write_template_size(record.template_size())?;
        } else if let Some(distance_to_next_fragment) = record.distance_to_next_fragment() {
            self.write_distance_to_next_fragment(distance_to_next_fragment)?;
        }

        Ok(())
    }

    fn write_next_mate_bit_flags(&mut self, next_mate_flags: NextMateFlags) -> io::Result<()> {
        let next_mate_bit_flags = i32::from(u8::from(next_mate_flags));

        self.compression_header
            .data_series_encoding_map()
            .next_mate_bit_flags_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::NextMateBitFlags),
                )
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                next_mate_bit_flags,
            )
    }

    fn write_next_fragment_reference_sequence_id(
        &mut self,
        next_fragment_reference_sequence_id: Option<usize>,
    ) -> io::Result<()> {
        const UNMAPPED: i32 = -1;

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .next_fragment_reference_sequence_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(
                        DataSeries::NextFragmentReferenceSequenceId,
                    ),
                )
            })?;

        let raw_next_fragment_reference_sequence_id =
            if let Some(id) = next_fragment_reference_sequence_id {
                i32::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?
            } else {
                UNMAPPED
            };

        encoding.encode(
            self.core_data_writer,
            self.external_data_writers,
            raw_next_fragment_reference_sequence_id,
        )
    }

    fn write_next_mate_alignment_start(
        &mut self,
        next_mate_alignment_start: Option<Position>,
    ) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .next_mate_alignment_start_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::NextMateAlignmentStart),
                )
            })?;

        let position = i32::try_from(
            next_mate_alignment_start
                .map(usize::from)
                .unwrap_or_default(),
        )
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, position)
    }

    fn write_template_size(&mut self, template_size: i32) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .template_size_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::TemplateSize),
                )
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                template_size,
            )
    }

    fn write_distance_to_next_fragment(
        &mut self,
        distance_to_next_fragment: usize,
    ) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .distance_to_next_fragment_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::DistanceToNextFragment),
                )
            })?;

        let n = i32::try_from(distance_to_next_fragment)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_tag_data(&mut self, record: &Record) -> io::Result<()> {
        use bam::record::codec::encoder::data::field::put_value;

        let preservation_map = self.compression_header.preservation_map();
        let tag_ids_dictionary = preservation_map.tag_ids_dictionary();

        let keys: Vec<_> = record
            .tags()
            .iter()
            .map(|(tag, value)| tag_ids_dictionary::Key::new(tag, value.ty()))
            .collect();

        let tag_line = tag_ids_dictionary
            .iter()
            .enumerate()
            .find(|(_, k)| **k == keys)
            .map(|(i, _)| i)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "tag line not in tag IDs dictionary",
                )
            })?;

        self.write_tag_line(tag_line)?;

        let tag_encoding_map = self.compression_header.tag_encoding_map();
        let mut buf = Vec::new();

        for result in sam::alignment::Record::data(record).iter() {
            let (tag, value) = result?;

            let key = tag_ids_dictionary::Key::new(tag, value.ty());
            let id = block::ContentId::from(key);
            let encoding = tag_encoding_map.get(&id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    WriteRecordError::MissingTagEncoding(key),
                )
            })?;

            buf.clear();
            put_value(&mut buf, &value)?;

            encoding.encode(self.core_data_writer, self.external_data_writers, &buf)?;
        }

        Ok(())
    }

    fn write_tag_line(&mut self, tag_line: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .tag_ids_encoding();

        let n =
            i32::try_from(tag_line).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mapped_read(&mut self, record: &Record) -> io::Result<()> {
        self.write_number_of_read_features(record.features().len())?;

        let mut prev_position = 0;

        for feature in record.features().iter() {
            let position = usize::from(feature.position()) - prev_position;
            self.write_feature(feature, position)?;
            prev_position = usize::from(feature.position());
        }

        self.write_mapping_quality(record.mapping_quality())?;

        if record.cram_flags().are_quality_scores_stored_as_array() {
            for &score in record.quality_scores().as_ref() {
                self.write_quality_score(score)?;
            }
        }

        Ok(())
    }

    fn write_number_of_read_features(&mut self, feature_count: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .number_of_read_features_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::NumberOfReadFeatures),
                )
            })?;

        let number_of_read_features = i32::try_from(feature_count)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(
            self.core_data_writer,
            self.external_data_writers,
            number_of_read_features,
        )
    }

    fn write_feature(&mut self, feature: &Feature, position: usize) -> io::Result<()> {
        self.write_feature_code(feature.code())?;
        self.write_feature_position(position)?;

        match feature {
            Feature::Bases(_, bases) => {
                self.write_stretches_of_bases(bases)?;
            }
            Feature::Scores(_, quality_scores) => {
                self.write_stretches_of_quality_scores(quality_scores)?;
            }
            Feature::ReadBase(_, base, quality_score) => {
                self.write_base(*base)?;
                self.write_quality_score(*quality_score)?;
            }
            Feature::Substitution(_, value) => {
                self.write_base_substitution_code(*value)?;
            }
            Feature::Insertion(_, bases) => {
                self.write_insertion(bases)?;
            }
            Feature::Deletion(_, len) => {
                self.write_deletion_length(*len)?;
            }
            Feature::InsertBase(_, base) => {
                self.write_base(*base)?;
            }
            Feature::QualityScore(_, score) => {
                self.write_quality_score(*score)?;
            }
            Feature::ReferenceSkip(_, len) => {
                self.write_reference_skip_length(*len)?;
            }
            Feature::SoftClip(_, bases) => {
                self.write_soft_clip(bases)?;
            }
            Feature::Padding(_, len) => {
                self.write_padding(*len)?;
            }
            Feature::HardClip(_, len) => {
                self.write_hard_clip(*len)?;
            }
        }

        Ok(())
    }

    fn write_feature_code(&mut self, code: feature::Code) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .read_features_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::ReadFeaturesCodes),
                )
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                u8::from(code),
            )
    }

    fn write_feature_position(&mut self, position: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .in_read_positions_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::InReadPositions),
                )
            })?;

        let position =
            i32::try_from(position).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, position)
    }

    fn write_stretches_of_bases(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .stretches_of_bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::StretchesOfBases),
                )
            })?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_stretches_of_quality_scores(&mut self, quality_scores: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .stretches_of_quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(
                        DataSeries::StretchesOfQualityScores,
                    ),
                )
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                quality_scores,
            )
    }

    fn write_base(&mut self, base: u8) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::Bases),
                )
            })?
            .encode(self.core_data_writer, self.external_data_writers, base)
    }

    fn write_quality_score(&mut self, quality_score: u8) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::QualityScores),
                )
            })?
            .encode(
                self.core_data_writer,
                self.external_data_writers,
                quality_score,
            )
    }

    fn write_base_substitution_code(&mut self, value: substitution::Value) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .base_substitution_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::BaseSubstitutionCodes),
                )
            })?;

        let substitution_matrix = self
            .compression_header
            .preservation_map()
            .substitution_matrix();

        let code = match value {
            substitution::Value::Bases(reference_base, read_base) => {
                substitution_matrix.find_code(reference_base, read_base)
            }
            substitution::Value::Code(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "base substitution cannot be a code on write",
                ));
            }
        };

        encoding.encode(self.core_data_writer, self.external_data_writers, code)
    }

    fn write_insertion(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .insertion_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::Insertion),
                )
            })?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_deletion_length(&mut self, len: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .deletion_lengths_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::DeletionLengths),
                )
            })?;

        let n = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_reference_skip_length(&mut self, len: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .reference_skip_length_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSkipLength),
                )
            })?;

        let n = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_soft_clip(&mut self, bases: &[u8]) -> io::Result<()> {
        self.compression_header
            .data_series_encoding_map()
            .soft_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::SoftClip),
                )
            })?
            .encode(self.core_data_writer, self.external_data_writers, bases)
    }

    fn write_padding(&mut self, len: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .padding_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::Padding),
                )
            })?;

        let n = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_hard_clip(&mut self, len: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .hard_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::HardClip),
                )
            })?;

        let n = i32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_mapping_quality(
        &mut self,
        mapping_quality: Option<sam::alignment::record_buf::MappingQuality>,
    ) -> io::Result<()> {
        const MISSING: u8 = 0xff;

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .mapping_qualities_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    WriteRecordError::MissingDataSeriesEncoding(DataSeries::MappingQualities),
                )
            })?;

        let n = i32::from(mapping_quality.map(u8::from).unwrap_or(MISSING));
        encoding.encode(self.core_data_writer, self.external_data_writers, n)
    }

    fn write_unmapped_read(&mut self, record: &Record) -> io::Result<()> {
        for &base in record.bases().as_ref() {
            self.write_base(base)?;
        }

        if record.cram_flags().are_quality_scores_stored_as_array() {
            for &score in record.quality_scores().as_ref() {
                self.write_quality_score(score)?;
            }
        }

        Ok(())
    }
}
