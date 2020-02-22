use std::io::{self, BufRead};

use byteorder::ReadBytesExt;

use crate::{
    data_series::DataSeries,
    encoding::{self, Encoding},
    huffman::CanonicalHuffmanDecoder,
    num::{read_itf8, Itf8},
    BitReader, Block, CompressionHeader, Record,
};

pub struct Reader<'a> {
    pub compression_header: &'a CompressionHeader,
    core_data_block: &'a Block,
    pub external_blocks: &'a [Block],
    reference_sequence_id: Itf8,
    prev_alignment_start: Itf8,
}

impl<'a> Reader<'a> {
    pub fn new(
        compression_header: &'a CompressionHeader,
        core_data_block: &'a Block,
        external_blocks: &'a [Block],
        reference_sequence_id: Itf8,
        initial_alignment_start: Itf8,
    ) -> Self {
        Self {
            compression_header,
            core_data_block,
            external_blocks,
            reference_sequence_id,
            prev_alignment_start: initial_alignment_start,
        }
    }

    pub fn read_record(&self, record: &mut Record) -> io::Result<()> {
        record.bam_bit_flags = self.read_bam_bit_flags()?;
        record.cram_bit_flags = self.read_cram_bit_flags()?;

        self.read_positional_data(record)?;

        let preservation_map = self.compression_header.preservation_map();

        if preservation_map.read_names_included() {
            record.read_name = self.read_read_name()?;
        } else {
            todo!("generate name");
        }

        Ok(())
    }

    fn read_bam_bit_flags(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::BamBitFlags)
            .expect("missing BF");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
    }

    fn read_cram_bit_flags(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::CramBitFlags)
            .expect("missing CF");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
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

        record.read_group = self.read_read_group()?;

        Ok(())
    }

    fn read_reference_id(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReferenceId)
            .expect("missing RI");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
    }

    fn read_read_length(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadLengths)
            .expect("missing RL");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
    }

    fn read_alignment_start(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::InSeqPositions)
            .expect("missing AP");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
    }

    fn read_read_group(&self) -> io::Result<Itf8> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadGroups)
            .expect("missing RG");

        decode_itf8(&encoding, &self.core_data_block, &self.external_blocks)
    }

    pub fn read_read_name(&self) -> io::Result<Vec<u8>> {
        let encoding = self
            .compression_header
            .data_series_encoding_map()
            .get(&DataSeries::ReadNames)
            .expect("missing RN");

        decode_byte_array(&encoding, &self.external_blocks)
    }
}

fn decode_itf8(
    encoding: &Encoding,
    core_data_block: &Block,
    external_blocks: &[Block],
) -> io::Result<Itf8> {
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
        encoding::Kind::Huffman => {
            let mut reader = encoding.args();

            let alphabet_len = read_itf8(&mut reader)? as usize;
            let mut alphabet = Vec::with_capacity(alphabet_len);

            for _ in 0..alphabet_len {
                let symbol = read_itf8(&mut reader)?;
                alphabet.push(symbol);
            }

            let bit_lens_len = read_itf8(&mut reader)? as usize;
            let mut bit_lens = Vec::with_capacity(bit_lens_len);

            for _ in 0..bit_lens_len {
                let len = read_itf8(&mut reader)?;
                bit_lens.push(len);
            }

            let decoder = CanonicalHuffmanDecoder::new(&alphabet, &bit_lens);
            let data = core_data_block.decompressed_data();
            let mut reader = BitReader::new(&data[..]);
            decoder.read(&mut reader)
        }
        _ => todo!(),
    }
}

fn decode_byte_array(encoding: &Encoding, external_blocks: &[Block]) -> io::Result<Vec<u8>> {
    match encoding.kind() {
        encoding::Kind::ByteArrayStop => {
            let mut reader = encoding.args();
            let stop_byte = reader.read_u8()?;
            let block_content_id = read_itf8(&mut reader)?;

            let block = external_blocks
                .iter()
                .find(|b| b.content_id() == block_content_id)
                .expect("could not find block");

            let data = block.decompressed_data();
            let mut reader = &data[..];
            let mut buf = Vec::new();
            reader.read_until(stop_byte, &mut buf)?;

            Ok(buf)
        }
        _ => todo!(),
    }
}
