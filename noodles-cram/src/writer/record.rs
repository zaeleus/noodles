mod tag;

use std::{
    collections::HashMap,
    io::{self, Write},
};

use byteorder::WriteBytesExt;

use noodles_sam as sam;

use crate::{
    container::{
        compression_header::{data_series_encoding_map::DataSeries, Encoding},
        CompressionHeader, ReferenceSequenceId,
    },
    num::{write_itf8, Itf8},
    record::{feature, Feature, Flags, NextMateFlags},
    BitWriter, Record,
};

pub struct Writer<'a, W, X> {
    compression_header: &'a CompressionHeader,
    core_data_writer: &'a mut BitWriter<W>,
    external_data_writers: &'a mut HashMap<Itf8, X>,
    reference_sequence_id: ReferenceSequenceId,
    prev_alignment_start: Itf8,
}

impl<'a, W, X> Writer<'a, W, X>
where
    W: Write,
    X: Write,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_writer: &'a mut BitWriter<W>,
        external_data_writers: &'a mut HashMap<Itf8, X>,
        reference_sequence_id: ReferenceSequenceId,
        initial_alignment_start: Itf8,
    ) -> Self {
        Self {
            compression_header,
            core_data_writer,
            external_data_writers,
            reference_sequence_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write_bam_bit_flags(record.bam_flags())?;
        self.write_cram_bit_flags(record.flags())?;

        self.write_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            self.write_read_name(record.read_name())?;
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

    fn write_bam_bit_flags(&mut self, bam_flags: sam::record::Flags) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::BamBitFlags)
            .expect("missing BF");

        let bam_bit_flags = i32::from(u16::from(bam_flags));

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            bam_bit_flags,
        )
    }

    fn write_cram_bit_flags(&mut self, flags: Flags) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::CramBitFlags)
            .expect("missing CF");

        let cram_bit_flags = i32::from(u8::from(flags));

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            cram_bit_flags,
        )
    }

    fn write_positional_data(&mut self, record: &Record) -> io::Result<()> {
        if self.reference_sequence_id.is_many() {
            let reference_id = i32::from(record.reference_sequence_id());
            self.write_reference_id(reference_id)?;
        }

        self.write_read_length(record.read_length())?;

        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let alignment_start = if ap_data_series_delta {
            record.alignment_start() - self.prev_alignment_start
        } else {
            record.alignment_start()
        };

        self.write_alignment_start(alignment_start)?;

        let read_group = i32::from(record.read_group_id());
        self.write_read_group(read_group)?;

        Ok(())
    }

    fn write_reference_id(&mut self, reference_id: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReferenceId)
            .expect("missing RI");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            reference_id,
        )
    }

    fn write_read_length(&mut self, read_length: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadLengths)
            .expect("missing RL");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            read_length,
        )
    }

    fn write_alignment_start(&mut self, alignment_start: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::InSeqPositions)
            .expect("missing AP");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            alignment_start,
        )
    }

    fn write_read_group(&mut self, read_group: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadGroups)
            .expect("missing RG");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            read_group,
        )
    }

    fn write_read_name(&mut self, read_name: &[u8]) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadNames)
            .expect("missing RN");

        encode_byte_array(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            read_name,
        )
    }

    fn write_mate_data(&mut self, record: &Record) -> io::Result<()> {
        let flags = record.flags();

        if flags.is_detached() {
            self.write_next_mate_bit_flags(record.next_mate_flags())?;

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                self.write_read_name(record.read_name())?;
            }

            let next_fragment_reference_sequence_id =
                i32::from(record.next_fragment_reference_sequence_id());
            self.write_next_fragment_reference_sequence_id(next_fragment_reference_sequence_id)?;

            self.write_next_mate_alignment_start(record.next_mate_alignment_start())?;
            self.write_template_size(record.template_size())?;
        } else {
            self.write_distance_to_next_fragment(record.distance_to_next_fragment())?;
        }

        Ok(())
    }

    fn write_next_mate_bit_flags(&mut self, next_mate_flags: NextMateFlags) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextMateBitFlags)
            .expect("missing MF");

        let next_mate_bit_flags = i32::from(u8::from(next_mate_flags));

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            next_mate_bit_flags,
        )
    }

    fn write_next_fragment_reference_sequence_id(
        &mut self,
        next_fragment_reference_sequence_id: Itf8,
    ) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextFragmentReferenceSequenceId)
            .expect("missing NS");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            next_fragment_reference_sequence_id,
        )
    }

    fn write_next_mate_alignment_start(
        &mut self,
        next_mate_alignment_start: Itf8,
    ) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NextMateAlignmentStart)
            .expect("missing NP");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            next_mate_alignment_start,
        )
    }

    fn write_template_size(&mut self, template_size: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::TemplateSize)
            .expect("missing TS");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            template_size,
        )
    }

    fn write_distance_to_next_fragment(
        &mut self,
        distance_to_next_fragment: Itf8,
    ) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::DistanceToNextFragment)
            .expect("missing NF");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            distance_to_next_fragment,
        )
    }

    fn write_tag_data(&mut self, record: &Record) -> io::Result<()> {
        let preservation_map = self.compression_header.preservation_map();
        let tag_ids_dictionary = preservation_map.tag_ids_dictionary();

        let keys: Vec<_> = record.tags().iter().map(|tag| tag.key()).collect();
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

        self.write_tag_line(tag_line as Itf8)?;

        let tag_encoding_map = self.compression_header.tag_encoding_map();

        for tag in record.tags() {
            let id = tag.key().id();
            let encoding = tag_encoding_map.get(&id).expect("missing tag encoding");

            let mut buf = Vec::new();
            tag::write_value(&mut buf, tag.value())?;

            encode_byte_array(
                encoding,
                &mut self.core_data_writer,
                &mut self.external_data_writers,
                &buf,
            )?;
        }

        Ok(())
    }

    fn write_tag_line(&mut self, tag_line: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::TagIds)
            .expect("missing TL");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            tag_line,
        )
    }

    fn write_mapped_read(&mut self, record: &Record) -> io::Result<()> {
        self.write_number_of_read_features(record.features().len())?;

        let mut prev_position = 0;

        for feature in record.features() {
            let position = feature.position() - prev_position;
            self.write_feature(feature, position)?;
            prev_position = feature.position();
        }

        let mapping_quality = i32::from(u8::from(record.mapping_quality()));
        self.write_mapping_quality(mapping_quality)?;

        let flags = record.flags();

        if flags.are_quality_scores_stored_as_array() {
            for &score in record.quality_scores() {
                self.write_quality_score(score)?;
            }
        }

        Ok(())
    }

    fn write_number_of_read_features(&mut self, feature_count: usize) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::NumberOfReadFeatures)
            .expect("missing FN");

        let number_of_read_features = feature_count as Itf8;

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            number_of_read_features,
        )
    }

    fn write_feature(&mut self, feature: &Feature, position: Itf8) -> io::Result<()> {
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
            Feature::Substitution(_, code) => {
                self.write_base_substitution_code(*code)?;
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
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadFeaturesCodes)
            .expect("missing FC");

        let feature_code = char::from(code) as Itf8;

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            feature_code,
        )
    }

    fn write_feature_position(&mut self, position: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::InReadPositions)
            .expect("missing FP");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            position,
        )
    }

    fn write_stretches_of_bases(&mut self, bases: &[u8]) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::StretchesOfBases)
            .expect("missing BB");

        encode_byte_array(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            bases,
        )
    }

    fn write_stretches_of_quality_scores(&mut self, quality_scores: &[u8]) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::StretchesOfQualityScores)
            .expect("missing QQ");

        encode_byte_array(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            quality_scores,
        )
    }

    fn write_base(&mut self, base: u8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::Bases)
            .expect("missing BA");

        encode_byte(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            base,
        )
    }

    fn write_quality_score(&mut self, quality_score: u8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::QualityScores)
            .expect("missing QS");

        encode_byte(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            quality_score,
        )
    }

    fn write_base_substitution_code(&mut self, code: u8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::BaseSubstitutionCodes)
            .expect("missing BS");

        encode_byte(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            code,
        )
    }

    fn write_insertion(&mut self, bases: &[u8]) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::Insertion)
            .expect("missing IN");

        encode_byte_array(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            bases,
        )
    }

    fn write_deletion_length(&mut self, len: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::DeletionLengths)
            .expect("missing DL");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            len,
        )
    }

    fn write_reference_skip_length(&mut self, len: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReferenceSkipLength)
            .expect("missing RS");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            len,
        )
    }

    fn write_soft_clip(&mut self, bases: &[u8]) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::SoftClip)
            .expect("missing SC");

        encode_byte_array(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            bases,
        )
    }

    fn write_padding(&mut self, len: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::Padding)
            .expect("missing PD");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            len,
        )
    }

    fn write_hard_clip(&mut self, len: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::HardClip)
            .expect("missing HC");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            len,
        )
    }

    fn write_mapping_quality(&mut self, mapping_quality: Itf8) -> io::Result<()> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::MappingQualities)
            .expect("missing MQ");

        encode_itf8(
            &encoding,
            &mut self.core_data_writer,
            &mut self.external_data_writers,
            mapping_quality,
        )
    }

    fn write_unmapped_read(&mut self, record: &Record) -> io::Result<()> {
        for &base in record.bases() {
            self.write_base(base)?;
        }

        let flags = record.flags();

        if flags.are_quality_scores_stored_as_array() {
            for &score in record.quality_scores() {
                self.write_quality_score(score)?;
            }
        }

        Ok(())
    }
}

fn encode_byte<W, X>(
    encoding: &Encoding,
    _core_data_writer: &mut BitWriter<W>,
    external_data_writers: &mut HashMap<Itf8, X>,
    value: u8,
) -> io::Result<()>
where
    W: Write,
    X: Write,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let writer = external_data_writers
                .get_mut(&block_content_id)
                .expect("could not find block");

            writer.write_u8(value)
        }
        _ => todo!("encode_byte: {:?}", encoding),
    }
}

fn encode_itf8<W, X>(
    encoding: &Encoding,
    _core_data_writer: &mut BitWriter<W>,
    external_data_writers: &mut HashMap<Itf8, X>,
    value: Itf8,
) -> io::Result<()>
where
    W: Write,
    X: Write,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let writer = external_data_writers
                .get_mut(&block_content_id)
                .expect("could not find block");

            write_itf8(writer, value)
        }
        _ => todo!("encode_itf8: {:?}", encoding),
    }
}

fn encode_byte_array<W, X>(
    encoding: &Encoding,
    core_data_writer: &mut BitWriter<W>,
    external_data_writers: &mut HashMap<Itf8, X>,
    data: &[u8],
) -> io::Result<()>
where
    W: Write,
    X: Write,
{
    match encoding {
        Encoding::External(block_content_id) => {
            let writer = external_data_writers
                .get_mut(&block_content_id)
                .expect("could not find block");

            writer.write_all(data)
        }
        Encoding::ByteArrayLen(len_encoding, value_encoding) => {
            let len = data.len() as Itf8;
            encode_itf8(&len_encoding, core_data_writer, external_data_writers, len)?;

            encode_byte_array(
                &value_encoding,
                core_data_writer,
                external_data_writers,
                data,
            )
        }
        Encoding::ByteArrayStop(stop_byte, block_content_id) => {
            let writer = external_data_writers
                .get_mut(&block_content_id)
                .expect("could not find block");

            writer.write_all(data)?;
            writer.write_u8(*stop_byte)?;

            Ok(())
        }
        _ => todo!("encode_byte_array: {:?}", encoding),
    }
}
