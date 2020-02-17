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
    reference_sequence_id: Itf8,
    prev_alignment_start: Itf8,
}

impl<'a> Reader<'a> {
    pub fn new(
        compression_header: &'a CompressionHeader,
        external_blocks: &'a [Block],
        reference_sequence_id: Itf8,
        initial_alignment_start: Itf8,
    ) -> Self {
        Self {
            compression_header,
            external_blocks,
            reference_sequence_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn read_record(&self, record: &mut Record) -> io::Result<()> {
        record.bam_bit_flags = self.read_bam_bit_flags()?;
        record.cram_bit_flags = self.read_cram_bit_flags()?;

        self.read_positional_data(record)?;

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

    fn read_positional_data(&self, record: &mut Record) -> io::Result<()> {
        record.reference_id = if self.reference_sequence_id == -2 {
            self.read_reference_id()?
        } else {
            self.reference_sequence_id
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

        Ok(())
    }

    fn read_reference_id(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReferenceId)
            .expect("missing RI");

        decode_itf8(&encoding, &self.external_blocks)
    }

    fn read_read_length(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadLengths)
            .expect("missing RL");

        decode_itf8(&encoding, &self.external_blocks)
    }

    fn read_alignment_start(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::InSeqPositions)
            .expect("missing AP");

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
