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
    record::Flags,
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
        self.write_bam_bit_flags(record.bam_bit_flags())?;
        self.write_cram_bit_flags(record.cram_bit_flags())?;

        self.write_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            self.write_read_name(&record.read_name)?;
        }

        self.prev_alignment_start = record.alignment_start();

        todo!();
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
            self.write_reference_id(record.reference_id)?;
        }

        self.write_read_length(record.read_length)?;

        let ap_data_series_delta = self
            .compression_header
            .preservation_map()
            .ap_data_series_delta();

        let alignment_start = if ap_data_series_delta {
            record.alignment_start - self.prev_alignment_start
        } else {
            record.alignment_start
        };

        self.write_alignment_start(alignment_start)?;

        self.write_read_group(record.read_group)?;

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
    _core_data_writer: &mut BitWriter<W>,
    external_data_writers: &mut HashMap<Itf8, X>,
    data: &[u8],
) -> io::Result<()>
where
    W: Write,
    X: Write,
{
    match encoding {
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
