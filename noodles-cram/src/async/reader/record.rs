use std::{
    collections::HashMap,
    convert::TryFrom,
    io::{self, Read},
};

use noodles_bam as bam;
use noodles_sam as sam;
use tokio::io::{AsyncBufRead, AsyncBufReadExt, AsyncRead};

use crate::{
    container::ReferenceSequenceId,
    data_container::{
        compression_header::{data_series_encoding_map::DataSeries, Encoding},
        CompressionHeader,
    },
    num::Itf8,
    r#async::reader::num::read_itf8,
    reader::record::ReadRecordError,
    record::{Builder, Flags, ReadGroupId},
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

        builder = self.read_positional_data(builder).await?;
        builder = self.read_read_names(builder).await?;

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

    async fn read_positional_data(&mut self, mut builder: Builder) -> io::Result<Builder> {
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

        Ok(builder)
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

    async fn read_read_length(&mut self) -> io::Result<Itf8> {
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
                    io::Error::new(io::ErrorKind::InvalidData, "missing external block")
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
