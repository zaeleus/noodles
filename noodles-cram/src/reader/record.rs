use std::io;

use crate::{
    data_series::DataSeries,
    encoding::{self, Encoding},
    num::{read_itf8, Itf8},
    Block, CompressionHeader, Record,
};

pub struct Reader<'a> {
    pub compression_header: &'a CompressionHeader,
    pub external_blocks: &'a [Block],
}

impl<'a> Reader<'a> {
    pub fn new(compression_header: &'a CompressionHeader, external_blocks: &'a [Block]) -> Self {
        Self {
            compression_header,
            external_blocks,
        }
    }

    pub fn read_record(&self, record: &mut Record) -> io::Result<()> {
        record.bam_bit_flags = self.read_bam_bit_flags()?;
        record.cram_bit_flags = self.read_cram_bit_flags()?;
        Ok(())
    }

    fn read_bam_bit_flags(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::BamBitFlags)
            .expect("missing BF");

        decode_itf8(&encoding, &self.external_blocks)
    }

    fn read_cram_bit_flags(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::CramBitFlags)
            .expect("missing CF");

        decode_itf8(&encoding, &self.external_blocks)
    }
}

fn decode_itf8(encoding: &Encoding, external_blocks: &[Block]) -> io::Result<Itf8> {
    match encoding.kind() {
        encoding::Kind::External => {
            let mut reader = encoding.args();
            let block_content_id = read_itf8(&mut reader)?;

            let block = external_blocks
                .iter()
                .find(|b| b.content_id() == block_content_id)
                .expect("could not find block");

            let data = block.decompressed_data();
            let mut reader = &data[..];
            read_itf8(&mut reader)
        }
        _ => todo!(),
    }
}
