use std::{collections::HashMap, io};

use byteorder::WriteBytesExt;
use bytes::Buf;

use crate::{
    container::block,
    data_container::compression_header::encoding::{Decode, Encode},
    huffman::CanonicalHuffmanDecoder,
    io::{BitReader, BitWriter},
    reader::record::ExternalDataReaders,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Byte {
    // block_content_id
    External(block::ContentId),
    // alphabet, bit_lens
    Huffman(Vec<i32>, Vec<u32>),
}

impl Decode for Byte {
    type Value = u8;

    fn decode<R, S>(
        &self,
        core_data_reader: &mut BitReader<R>,
        external_data_readers: &mut ExternalDataReaders<S>,
    ) -> io::Result<Self::Value>
    where
        R: Buf,
        S: Buf,
    {
        match self {
            Byte::External(block_content_id) => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                if !src.has_remaining() {
                    return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
                }

                Ok(src.get_u8())
            }
            Byte::Huffman(alphabet, bit_lens) => {
                if alphabet.len() == 1 {
                    Ok(alphabet[0] as u8)
                } else {
                    let decoder = CanonicalHuffmanDecoder::new(alphabet, bit_lens);
                    decoder.decode(core_data_reader).map(|i| i as u8)
                }
            }
        }
    }
}

impl Encode for Byte {
    type Value = u8;

    fn encode<W, X>(
        &self,
        _core_data_writer: &mut BitWriter<W>,
        external_data_writers: &mut HashMap<block::ContentId, X>,
        value: &Self::Value,
    ) -> io::Result<()>
    where
        W: io::Write,
        X: io::Write,
    {
        match self {
            Byte::External(block_content_id) => {
                let writer = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                writer.write_u8(*value)
            }
            _ => todo!("encode_byte: {:?}", self),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_container::compression_header::Encoding;

    #[test]
    fn test_decode() -> io::Result<()> {
        fn t(encoding: &Encoding<Byte>, expected: u8) -> io::Result<()> {
            let core_data = [0b10000000];
            let mut core_data_reader = BitReader::new(&core_data[..]);

            let external_data = [0x0d];
            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(block::ContentId::from(1), &external_data[..]);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        t(
            &Encoding::new(Byte::External(block::ContentId::from(1))),
            0x0d,
        )?;
        t(&Encoding::new(Byte::Huffman(vec![0x4e], vec![0])), 0x4e)?;

        Ok(())
    }
}
