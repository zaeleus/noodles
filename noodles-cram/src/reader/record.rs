use std::io;

use crate::{
    data_series::DataSeries,
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
        record.bam_bit_flags = dbg!(self.read_bam_bit_flags())?;
        Ok(())
    }

    fn read_bam_bit_flags(&self) -> io::Result<Itf8> {
        let data_series_encoding_map = self.compression_header.data_series_encoding_map();

        let encoding = data_series_encoding_map
            .get(&DataSeries::BamBitFlags)
            .expect("missing BF");

        let mut reader = encoding.args();
        let block_content_id = dbg!(read_itf8(&mut reader))?;

        let block = self
            .external_blocks
            .iter()
            .find(|b| b.content_id() == block_content_id)
            .expect("could not find block");

        let data = block.decompressed_data();
        let mut reader = &data[..];
        read_itf8(&mut reader)
    }
}
