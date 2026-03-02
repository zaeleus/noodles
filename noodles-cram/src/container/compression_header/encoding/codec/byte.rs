use std::{borrow::Cow, io};

use crate::{
    container::{
        block,
        compression_header::encoding::{Decode, Encode},
    },
    huffman::{CanonicalHuffmanDecoder, CanonicalHuffmanEncoder},
    io::{
        BitReader, BitWriter, reader::container::slice::records::ExternalDataReaders,
        writer::container::slice::records::ExternalDataWriters,
    },
};

#[derive(Clone, Debug)]
pub enum Byte {
    External {
        block_content_id: block::ContentId,
    },
    Huffman {
        alphabet: Vec<i32>,
        bit_lens: Vec<u32>,
        decoder: CanonicalHuffmanDecoder,
        encoder: CanonicalHuffmanEncoder,
    },
    // CRAM 4.0 codec
    Constant {
        value: u8,
    },
}

impl Byte {
    pub fn huffman(alphabet: Vec<i32>, bit_lens: Vec<u32>) -> Self {
        let decoder = CanonicalHuffmanDecoder::new(&alphabet, &bit_lens);
        let encoder = CanonicalHuffmanEncoder::new(&alphabet, &bit_lens);
        Self::Huffman {
            alphabet,
            bit_lens,
            decoder,
            encoder,
        }
    }
}

impl PartialEq for Byte {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (
                Self::External {
                    block_content_id: a,
                },
                Self::External {
                    block_content_id: b,
                },
            ) => a == b,
            (
                Self::Huffman {
                    alphabet: a1,
                    bit_lens: a2,
                    ..
                },
                Self::Huffman {
                    alphabet: b1,
                    bit_lens: b2,
                    ..
                },
            ) => a1 == b1 && a2 == b2,
            (Self::Constant { value: a }, Self::Constant { value: b }) => a == b,
            _ => false,
        }
    }
}

impl Eq for Byte {}

impl Byte {
    pub fn decode_take<'de>(
        &self,
        core_data_reader: &mut BitReader<'de>,
        external_data_readers: &mut ExternalDataReaders<'de>,
        len: usize,
    ) -> io::Result<Cow<'de, [u8]>> {
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

                let (buf, rest) = src
                    .split_at_checked(len)
                    .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

                *src = rest;

                Ok(Cow::Borrowed(buf))
            }
            Self::Huffman {
                alphabet, decoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(Cow::Owned(vec![alphabet[0] as u8; len]))
                } else {
                    let mut buf = Vec::with_capacity(len);

                    for _ in 0..len {
                        let value = decoder.decode(core_data_reader)?;
                        buf.push(value as u8);
                    }

                    Ok(Cow::Owned(buf))
                }
            }
            Self::Constant { value } => Ok(Cow::Owned(vec![*value; len])),
        }
    }

    pub fn encode_extend(
        &self,
        core_data_writer: &mut BitWriter,
        external_data_writers: &mut ExternalDataWriters,
        src: &[u8],
    ) -> io::Result<()> {
        match self {
            Self::External { block_content_id } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                dst.extend(src);

                Ok(())
            }
            Self::Huffman {
                alphabet, encoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(())
                } else {
                    for &b in src {
                        encoder.encode(core_data_writer, b as i32)?;
                    }

                    Ok(())
                }
            }
            Self::Constant { .. } => Ok(()),
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
            Self::Huffman {
                alphabet, decoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(alphabet[0] as u8)
                } else {
                    decoder.decode(core_data_reader).map(|i| i as u8)
                }
            }
            Self::Constant { value } => Ok(*value),
        }
    }
}

impl Encode<'_> for Byte {
    type Value = u8;

    fn encode(
        &self,
        core_data_writer: &mut BitWriter,
        external_data_writers: &mut ExternalDataWriters,
        value: Self::Value,
    ) -> io::Result<()> {
        match self {
            Self::External { block_content_id } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                dst.push(value);

                Ok(())
            }
            Self::Huffman {
                alphabet, encoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(())
                } else {
                    encoder.encode(core_data_writer, value as i32)
                }
            }
            Self::Constant { .. } => Ok(()),
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

        assert_eq!(&*dst, &external_data[..]);

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
        t(&Encoding::new(Byte::huffman(vec![0x4e], vec![0])), 0x4e)?;
        t(&Encoding::new(Byte::Constant { value: 0xff }), 0xff)?;

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
            let mut core_data_writer = BitWriter::default();

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
