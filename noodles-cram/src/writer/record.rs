use std::{
    collections::HashMap,
    io::{self, Write},
};

use noodles_sam as sam;

use crate::{
    container::{
        compression_header::{data_series_encoding_map::DataSeries, encoding, Encoding},
        CompressionHeader,
    },
    num::{read_itf8, write_itf8, Itf8},
    record::Flags,
    BitWriter, Record,
};

pub struct Writer<'a, W, X> {
    compression_header: &'a CompressionHeader,
    core_data_writer: &'a mut BitWriter<W>,
    external_data_writers: &'a mut HashMap<Itf8, X>,
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
    ) -> Self {
        Self {
            compression_header,
            core_data_writer,
            external_data_writers,
        }
    }

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write_bam_bit_flags(record.bam_bit_flags())?;
        self.write_cram_bit_flags(record.cram_bit_flags())?;

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
    match encoding.kind() {
        encoding::Kind::External => {
            let mut reader = encoding.args();
            let block_content_id = read_itf8(&mut reader)?;

            let writer = external_data_writers
                .get_mut(&block_content_id)
                .expect("could not find block");

            write_itf8(writer, value)
        }
        _ => todo!("encode_itf8: unhandled encoding {:?}", encoding.kind()),
    }
}
