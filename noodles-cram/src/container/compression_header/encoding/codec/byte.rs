use std::{collections::HashMap, io};

use byteorder::WriteBytesExt;

use crate::{
    container::{
        block,
        compression_header::encoding::{Decode, Encode},
    },
    huffman::CanonicalHuffmanDecoder,
    io::{
        reader::{container::slice::records::ExternalDataReaders, split_at_checked},
        BitReader, BitWriter,
    },
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Byte {
    External {
        block_content_id: block::ContentId,
    },
    Huffman {
        alphabet: Vec<i32>,
        bit_lens: Vec<u32>,
    },
}

impl Byte {
    pub fn decode_take<'de>(
        &self,
        _core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
        len: usize,
    ) -> io::Result<&'de [u8]> {
        match self {
            Self::External { block_content_id } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let (buf, rest) = split_at_checked(src, len)
                    .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

                *src = rest;

                Ok(buf)
            }
            Self::Huffman { .. } => todo!(),
        }
    }
}

impl<'de> Decode<'de> for Byte {
    type Value = u8;

    fn decode(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
    ) -> io::Result<Self::Value> {
        match self {
            Self::External { block_content_id } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let (n, rest) = src
                    .split_first()
                    .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

                *src = rest;

                Ok(*n)
            }
            Self::Huffman { alphabet, bit_lens } => {
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

impl Encode<'_> for Byte {
    type Value = u8;

    fn encode<W, X>(
        &self,
        _core_data_writer: &mut BitWriter<W>,
        external_data_writers: &mut HashMap<block::ContentId, X>,
        value: Self::Value,
    ) -> io::Result<()>
    where
        W: io::Write,
        X: io::Write,
    {
        match self {
            Self::External { block_content_id } => {
                let writer = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                writer.write_u8(value)
            }
            _ => todo!("encode_byte: {:?}", self),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::compression_header::Encoding;

    #[test]
    fn test_decode_take() -> io::Result<()> {
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let external_data = b"ndls";
        let mut external_data_readers = ExternalDataReaders::new();
        external_data_readers.insert(1, &external_data[..]);

        let codec = Byte::External {
            block_content_id: 1,
        };
        let dst = codec.decode_take(&mut core_data_reader, &mut external_data_readers, 4)?;

        assert_eq!(dst, external_data);

        Ok(())
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        fn t(encoding: &Encoding<Byte>, expected: u8) -> io::Result<()> {
            let core_data = [0b10000000];
            let mut core_data_reader = BitReader::new(&core_data[..]);

            let external_data = [0x0d];
            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(1, &external_data[..]);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        t(
            &Encoding::new(Byte::External {
                block_content_id: 1,
            }),
            0x0d,
        )?;
        t(
            &Encoding::new(Byte::Huffman {
                alphabet: vec![0x4e],
                bit_lens: vec![0],
            }),
            0x4e,
        )?;

        Ok(())
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        fn t(
            encoding: &Encoding<Byte>,
            value: u8,
            expected_core_data: &[u8],
            expected_external_data: &[u8],
        ) -> io::Result<()> {
            let mut core_data_writer = BitWriter::new(Vec::new());

            let block_content_id = 1;
            let mut external_data_writers = [(block_content_id, Vec::new())].into_iter().collect();

            encoding.encode(&mut core_data_writer, &mut external_data_writers, value)?;

            let actual_core_data = core_data_writer.finish()?;
            assert_eq!(actual_core_data, expected_core_data);

            let actual_external_data = &external_data_writers[&block_content_id];
            assert_eq!(actual_external_data, expected_external_data);

            Ok(())
        }

        t(
            &Encoding::new(Byte::External {
                block_content_id: 1,
            }),
            0x0d,
            &[],
            &[0x0d],
        )?;

        Ok(())
    }
}
