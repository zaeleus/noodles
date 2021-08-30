use std::{
    collections::HashMap,
    convert::TryFrom,
    io::{self, Read},
};

use noodles_bam as bam;
use noodles_sam as sam;
use tokio::io::{AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt};

use crate::{
    container::ReferenceSequenceId,
    data_container::{
        compression_header::{data_series_encoding_map::DataSeries, Encoding},
        CompressionHeader,
    },
    num::Itf8,
    r#async::reader::num::read_itf8,
    reader::record::ReadRecordError,
    record::{feature, Builder, Feature, Flags, NextMateFlags, ReadGroupId, Tag},
    BitReader, Record,
};

pub struct Reader<'a, CDR, EDR>
where
    CDR: Read,
    EDR: AsyncBufRead + Unpin,
{
    compression_header: &'a CompressionHeader,
    core_data_reader: BitReader<CDR>,
    external_data_readers: HashMap<Itf8, EDR>,
    reference_sequence_id: ReferenceSequenceId,
    prev_alignment_start: Itf8,
}

impl<'a, CDR, EDR> Reader<'a, CDR, EDR>
where
    CDR: Read,
    EDR: AsyncBufRead + Unpin,
{
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_reader: BitReader<CDR>,
        external_data_readers: HashMap<Itf8, EDR>,
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

    pub async fn read_record(&mut self) -> io::Result<Record> {
        let mut builder = Record::builder();

        let bam_bit_flags = self.read_bam_bit_flags().await?;
        builder = builder.set_bam_flags(bam_bit_flags);

        let cram_bit_flags = self.read_cram_bit_flags().await?;
        builder = builder.set_flags(cram_bit_flags);

        let (mut builder, read_length) = self.read_positional_data(builder).await?;
        builder = self.read_read_names(builder).await?;
        builder = self
            .read_mate_data(builder, bam_bit_flags, cram_bit_flags)
            .await?;

        let tags = self.read_tag_data().await?;
        builder = builder.set_tags(tags);

        builder = if bam_bit_flags.is_unmapped() {
            self.read_unmapped_read(builder, cram_bit_flags, read_length)
                .await?
        } else {
            self.read_mapped_read(builder, cram_bit_flags, read_length)
                .await?
        };

        let record = builder.build();

        self.prev_alignment_start = record.alignment_start();

        Ok(record)
    }

    async fn read_bam_bit_flags(&mut self) -> io::Result<sam::record::Flags> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .bam_bit_flags_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| u16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(sam::record::Flags::from)
    }

    async fn read_cram_bit_flags(&mut self) -> io::Result<Flags> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .cram_bit_flags_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(Flags::from)
    }

    async fn read_positional_data(&mut self, mut builder: Builder) -> io::Result<(Builder, usize)> {
        let reference_id = if self.reference_sequence_id.is_many() {
            self.read_reference_id().await?
        } else {
            i32::from(self.reference_sequence_id)
        };

        if reference_id != bam::record::reference_sequence_id::UNMAPPED {
            let reference_sequence_id = bam::record::ReferenceSequenceId::try_from(reference_id)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        let read_length = self.read_read_length().await?;
        builder = builder.set_read_length(read_length);

        let alignment_start = self.read_alignment_start().await?;
        builder = builder.set_alignment_start(alignment_start);

        let read_group = self.read_read_group().await?;
        builder = builder.set_read_group_id(read_group);

        Ok((builder, read_length))
    }

    async fn read_reference_id(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .reference_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceId),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_read_length(&mut self) -> io::Result<usize> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_lengths_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
    }

    async fn read_alignment_start(&mut self) -> io::Result<Itf8> {
        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .in_seq_positions_encoding();

        let alignment_start = decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await?;

        if ap_data_series_delta {
            Ok(self.prev_alignment_start + alignment_start)
        } else {
            Ok(alignment_start)
        }
    }

    async fn read_read_group(&mut self) -> io::Result<ReadGroupId> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_groups_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .map(ReadGroupId::from)
    }

    async fn read_read_names(&mut self, mut builder: Builder) -> io::Result<Builder> {
        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            let read_name = self.read_read_name().await?;
            builder = builder.set_read_name(read_name);
        }

        Ok(builder)
    }

    async fn read_read_name(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_names_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadNames),
                )
            })?;

        decode_byte_array(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_mate_data(
        &mut self,
        mut builder: Builder,
        mut bam_flags: sam::record::Flags,
        flags: Flags,
    ) -> io::Result<Builder> {
        if flags.is_detached() {
            let next_mate_bit_flags = self.read_next_mate_bit_flags().await?;
            builder = builder.set_next_mate_flags(next_mate_bit_flags);

            if next_mate_bit_flags.is_on_negative_strand() {
                bam_flags |= sam::record::Flags::MATE_REVERSE_COMPLEMENTED;
            }

            if next_mate_bit_flags.is_unmapped() {
                bam_flags |= sam::record::Flags::MATE_UNMAPPED;
            }

            builder = builder.set_bam_flags(bam_flags);

            let preservation_map = self.compression_header.preservation_map();

            if !preservation_map.read_names_included() {
                let read_name = self.read_read_name().await?;
                builder = builder.set_read_name(read_name);
            }

            if let Some(id) = self.read_next_fragment_reference_sequence_id().await? {
                builder = builder.set_next_fragment_reference_sequence_id(id);
            }

            let next_mate_alignment_start = self.read_next_mate_alignment_start().await?;
            builder = builder.set_next_mate_alignment_start(next_mate_alignment_start);

            let template_size = self.read_template_size().await?;
            builder = builder.set_template_size(template_size);
        } else if flags.has_mate_downstream() {
            let distance_to_next_fragment = self.read_distance_to_next_fragment().await?;
            builder = builder.set_distance_to_next_fragment(distance_to_next_fragment);
        }

        Ok(builder)
    }

    async fn read_next_mate_bit_flags(&mut self) -> io::Result<NextMateFlags> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .next_mate_bit_flags_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateBitFlags),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(NextMateFlags::from)
    }

    async fn read_next_fragment_reference_sequence_id(
        &mut self,
    ) -> io::Result<Option<bam::record::ReferenceSequenceId>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .next_fragment_reference_sequence_id_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(
                        DataSeries::NextFragmentReferenceSequenceId,
                    ),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|id| {
            if id == bam::record::reference_sequence_id::UNMAPPED {
                Ok(None)
            } else {
                bam::record::ReferenceSequenceId::try_from(id)
                    .map(Some)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }
        })
    }

    async fn read_next_mate_alignment_start(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .next_mate_alignment_start_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NextMateAlignmentStart),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_template_size(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .template_size_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::TemplateSize),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_distance_to_next_fragment(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .distance_to_next_fragment_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DistanceToNextFragment),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_tag_data(&mut self) -> io::Result<Vec<Tag>> {
        let tag_line = self.read_tag_line().await.and_then(|i| {
            usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let tag_keys = self
            .compression_header
            .preservation_map()
            .tag_ids_dictionary()
            .get(tag_line)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid tag line"))?;

        let tag_encoding_map = self.compression_header.tag_encoding_map();

        let mut tags = Vec::with_capacity(tag_keys.len());

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
            )
            .await?;

            let mut data_reader = bam::record::data::Reader::new(&data[..]);
            let value = data_reader.read_value_type(key.ty())?;

            let tag = Tag::new(*key, value);
            tags.push(tag);
        }

        Ok(tags)
    }

    async fn read_tag_line(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .tag_ids_encoding();

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_mapped_read(
        &mut self,
        mut builder: Builder,
        flags: Flags,
        read_length: usize,
    ) -> io::Result<Builder> {
        let feature_count = self.read_number_of_read_features().await?;

        let mut prev_position = 0;

        for _ in 0..feature_count {
            let feature = self.read_feature(prev_position).await?;
            prev_position = feature.position();
            builder = builder.add_feature(feature);
        }

        let mapping_quality = self.read_mapping_quality().await?;
        builder = builder.set_mapping_quality(mapping_quality);

        if flags.are_quality_scores_stored_as_array() {
            for _ in 0..read_length {
                let quality_score = self.read_quality_score().await?;
                builder = builder.add_quality_score(quality_score);
            }
        }

        Ok(builder)
    }

    async fn read_number_of_read_features(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .number_of_read_features_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::NumberOfReadFeatures),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_feature(&mut self, prev_position: i32) -> io::Result<Feature> {
        use feature::Code;

        let code = self.read_feature_code().await?;

        let delta = self.read_feature_position().await?;
        let position = prev_position + delta;

        match code {
            Code::Bases => {
                let bases = self.read_stretches_of_bases().await?;
                Ok(Feature::Bases(position, bases))
            }
            Code::Scores => {
                let quality_scores = self.read_stretches_of_quality_scores().await?;
                Ok(Feature::Scores(position, quality_scores))
            }
            Code::ReadBase => {
                let base = self.read_base().await?;
                let quality_score = self.read_quality_score().await?;
                Ok(Feature::ReadBase(position, base, quality_score))
            }
            Code::Substitution => {
                let code = self.read_base_substitution_code().await?;
                Ok(Feature::Substitution(position, code))
            }
            Code::Insertion => {
                let bases = self.read_insertion().await?;
                Ok(Feature::Insertion(position, bases))
            }
            Code::Deletion => {
                let len = self.read_deletion_length().await?;
                Ok(Feature::Deletion(position, len))
            }
            Code::InsertBase => {
                let base = self.read_base().await?;
                Ok(Feature::InsertBase(position, base))
            }
            Code::QualityScore => {
                let quality_score = self.read_quality_score().await?;
                Ok(Feature::QualityScore(position, quality_score))
            }
            Code::ReferenceSkip => {
                let len = self.read_reference_skip_length().await?;
                Ok(Feature::ReferenceSkip(position, len))
            }
            Code::SoftClip => {
                let bases = self.read_soft_clip().await?;
                Ok(Feature::SoftClip(position, bases))
            }
            Code::Padding => {
                let len = self.read_padding().await?;
                Ok(Feature::Padding(position, len))
            }
            Code::HardClip => {
                let len = self.read_hard_clip().await?;
                Ok(Feature::HardClip(position, len))
            }
        }
    }

    async fn read_feature_code(&mut self) -> io::Result<feature::Code> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .read_features_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReadFeaturesCodes),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .map(|id| id as u8 as char)
        .and_then(|id| {
            feature::Code::try_from(id).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    async fn read_feature_position(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .in_read_positions_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::InReadPositions),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_stretches_of_bases(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .stretches_of_bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::StretchesOfBases),
                )
            })?;

        decode_byte_array(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_stretches_of_quality_scores(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .stretches_of_quality_scores_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(
                        DataSeries::StretchesOfQualityScores,
                    ),
                )
            })?;

        decode_byte_array(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_base(&mut self) -> io::Result<u8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .bases_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Bases),
                )
            })?;

        decode_byte(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_quality_score(&mut self) -> io::Result<u8> {
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

        decode_byte(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_base_substitution_code(&mut self) -> io::Result<u8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .base_substitution_codes_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::BaseSubstitutionCodes),
                )
            })?;

        decode_byte(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_insertion(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .insertion_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Insertion),
                )
            })?;

        decode_byte_array(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_deletion_length(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .deletion_lengths_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::DeletionLengths),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_reference_skip_length(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .reference_skip_length_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::ReferenceSkipLength),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_soft_clip(&mut self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .soft_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::SoftClip),
                )
            })?;

        decode_byte_array(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_padding(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .padding_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::Padding),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_hard_clip(&mut self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .hard_clip_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::HardClip),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
    }

    async fn read_mapping_quality(&mut self) -> io::Result<sam::record::MappingQuality> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .mapping_qualities_encoding()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    ReadRecordError::MissingDataSeriesEncoding(DataSeries::MappingQualities),
                )
            })?;

        decode_itf8(
            encoding,
            &mut self.core_data_reader,
            &mut self.external_data_readers,
        )
        .await
        .and_then(|n| u8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .map(sam::record::MappingQuality::from)
    }

    async fn read_unmapped_read(
        &mut self,
        mut builder: Builder,
        flags: Flags,
        read_length: usize,
    ) -> io::Result<Builder> {
        for _ in 0..read_length {
            let base = self.read_base().await?;
            builder = builder.add_base(base);
        }

        if flags.are_quality_scores_stored_as_array() {
            for _ in 0..read_length {
                let quality_score = self.read_quality_score().await?;
                builder = builder.add_quality_score(quality_score);
            }
        }

        Ok(builder)
    }
}

async fn decode_byte<CDR, EDR>(
    encoding: &Encoding,
    _core_data_reader: &mut BitReader<CDR>,
    external_data_readers: &mut HashMap<Itf8, EDR>,
) -> io::Result<u8>
where
    CDR: Read,
    EDR: AsyncRead + Unpin,
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

            reader.read_u8().await
        }
        _ => todo!("decode_byte: {:?}", encoding),
    }
}

async fn decode_itf8<CDR, EDR>(
    encoding: &Encoding,
    _core_data_reader: &mut BitReader<CDR>,
    external_data_readers: &mut HashMap<Itf8, EDR>,
) -> io::Result<Itf8>
where
    CDR: Read,
    EDR: AsyncRead + Unpin,
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

            read_itf8(reader).await
        }
        _ => todo!("decode_itf8: {:?}", encoding),
    }
}

async fn decode_byte_array<CDR, EDR>(
    encoding: &Encoding,
    _core_data_reader: &mut BitReader<CDR>,
    external_data_readers: &mut HashMap<Itf8, EDR>,
) -> io::Result<Vec<u8>>
where
    CDR: Read,
    EDR: AsyncBufRead + Unpin,
{
    match encoding {
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
            reader.read_until(*stop_byte, &mut buf).await?;

            // Remove stop byte.
            buf.pop();

            Ok(buf)
        }
        _ => todo!("decode_byte_array: {:?}", encoding),
    }
}
