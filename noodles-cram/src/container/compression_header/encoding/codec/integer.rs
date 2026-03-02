use std::io;

use crate::{
    container::{
        block,
        compression_header::encoding::{Decode, Encode},
    },
    huffman::{CanonicalHuffmanDecoder, CanonicalHuffmanEncoder},
    io::{
        BitReader, BitWriter,
        reader::{
            container::slice::records::ExternalDataReaders,
            num::{read_itf8, read_sint7_64, read_uint7_64},
        },
        writer::{
            container::slice::records::ExternalDataWriters,
            num::{write_itf8, write_sint7_64, write_uint7_64},
        },
    },
};

#[derive(Clone, Debug)]
pub enum Integer {
    External {
        block_content_id: block::ContentId,
    },
    Golomb {
        offset: i32,
        m: i32,
    },
    Huffman {
        alphabet: Vec<i32>,
        bit_lens: Vec<u32>,
        decoder: CanonicalHuffmanDecoder,
        encoder: CanonicalHuffmanEncoder,
    },
    Beta {
        offset: i32,
        len: u32,
    },
    Subexp {
        offset: i32,
        k: i32,
    },
    GolombRice {
        offset: i32,
        log2_m: i32,
    },
    Gamma {
        offset: i32,
    },
    // CRAM 4.0 codecs
    VarintUnsigned {
        block_content_id: block::ContentId,
        offset: i64,
    },
    VarintSigned {
        block_content_id: block::ContentId,
        offset: i64,
    },
    ConstInt {
        value: i32,
    },
}

impl Integer {
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

impl PartialEq for Integer {
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
            (Self::Golomb { offset: a1, m: a2 }, Self::Golomb { offset: b1, m: b2 }) => {
                a1 == b1 && a2 == b2
            }
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
            (
                Self::Beta {
                    offset: a1,
                    len: a2,
                },
                Self::Beta {
                    offset: b1,
                    len: b2,
                },
            ) => a1 == b1 && a2 == b2,
            (Self::Subexp { offset: a1, k: a2 }, Self::Subexp { offset: b1, k: b2 }) => {
                a1 == b1 && a2 == b2
            }
            (
                Self::GolombRice {
                    offset: a1,
                    log2_m: a2,
                },
                Self::GolombRice {
                    offset: b1,
                    log2_m: b2,
                },
            ) => a1 == b1 && a2 == b2,
            (Self::Gamma { offset: a }, Self::Gamma { offset: b }) => a == b,
            (
                Self::VarintUnsigned {
                    block_content_id: a1,
                    offset: a2,
                },
                Self::VarintUnsigned {
                    block_content_id: b1,
                    offset: b2,
                },
            ) => a1 == b1 && a2 == b2,
            (
                Self::VarintSigned {
                    block_content_id: a1,
                    offset: a2,
                },
                Self::VarintSigned {
                    block_content_id: b1,
                    offset: b2,
                },
            ) => a1 == b1 && a2 == b2,
            (Self::ConstInt { value: a }, Self::ConstInt { value: b }) => a == b,
            _ => false,
        }
    }
}

impl Eq for Integer {}

impl<'de> Decode<'de> for Integer {
    type Value = i64;

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

                read_itf8(src).map(i64::from)
            }
            Self::Huffman {
                alphabet, decoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(i64::from(alphabet[0]))
                } else {
                    decoder.decode(core_data_reader).map(i64::from)
                }
            }
            Self::Beta { offset, len } => core_data_reader
                .read_i32(*len)
                .map(|i| i64::from(i - offset)),
            Self::Gamma { offset } => {
                let mut n = 0;

                while core_data_reader.read_bit()? == 0 {
                    n += 1;
                }

                let m = core_data_reader.read_i32(n)?;
                let x = (1 << n) + m;

                Ok(i64::from(x - offset))
            }
            Self::VarintUnsigned {
                block_content_id,
                offset,
            } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                read_uint7_64(src)
                    .and_then(|n| {
                        i64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    })
                    .map(|n| n + offset)
            }
            Self::VarintSigned {
                block_content_id,
                offset,
            } => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                read_sint7_64(src).map(|n| n + offset)
            }
            Self::ConstInt { value } => Ok(i64::from(*value)),
            Self::Golomb { offset, m } => {
                if *m <= 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid Golomb parameter: m={}", m),
                    ));
                }

                let mut q = 0i32;
                while core_data_reader.read_bit()? == 0 {
                    q += 1;
                }

                let b = 32 - ((*m - 1).leading_zeros());

                let value = if b == 0 {
                    q
                } else {
                    let r = core_data_reader.read_i32(b - 1)?;
                    let threshold = (1i32 << b) - m;

                    if r < threshold {
                        q * m + r
                    } else {
                        let r = (r << 1) | core_data_reader.read_i32(1)?;
                        q * m + r - threshold
                    }
                };

                Ok(i64::from(value - offset))
            }
            Self::Subexp { offset, k } => {
                if *k < 0 || *k >= 32 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid subexp parameter: k={k}"),
                    ));
                }
                let k = *k as u32;
                let mut n = 0u32;

                while core_data_reader.read_bit()? == 1 {
                    n += 1;
                }

                let value = if n < 2 {
                    core_data_reader.read_i32(k + n)?
                } else {
                    let extra = core_data_reader.read_i32(k + n + n - 1)?;
                    extra + (1 << (k + 1)) - (1 << k)
                };

                Ok(i64::from(value - offset))
            }
            Self::GolombRice { offset, log2_m } => {
                if *log2_m < 0 || *log2_m >= 32 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid Golomb-Rice parameter: log2_m={log2_m}"),
                    ));
                }
                let log2_m = *log2_m as u32;

                let mut q = 0i32;
                while core_data_reader.read_bit()? == 0 {
                    q += 1;
                }

                let r = core_data_reader.read_i32(log2_m)?;
                let value = (q << log2_m) | r;

                Ok(i64::from(value - offset))
            }
        }
    }
}

impl Encode<'_> for Integer {
    type Value = i64;

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

                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                write_itf8(dst, value)
            }
            Self::Huffman {
                alphabet, encoder, ..
            } => {
                if alphabet.len() == 1 {
                    Ok(())
                } else {
                    let value = i32::try_from(value)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                    encoder.encode(core_data_writer, value)
                }
            }
            Self::ConstInt { .. } => Ok(()),
            Self::Beta { offset, len } => {
                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                core_data_writer.write_u32((value + offset) as u32, *len as usize)
            }
            Self::Gamma { offset } => {
                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let x = value + offset;
                if x < 1 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("gamma encoding requires x >= 1, got {x}"),
                    ));
                }
                let n = 31 - (x as u32).leading_zeros();
                for _ in 0..n {
                    core_data_writer.write_u32(0, 1)?;
                }
                core_data_writer.write_u32(1, 1)?;
                let m = x - (1 << n);
                core_data_writer.write_u32(m as u32, n as usize)
            }
            Self::Golomb { offset, m } => {
                if *m <= 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid Golomb parameter: m={}", m),
                    ));
                }

                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let n = value + offset;

                if n < 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Golomb encoding requires non-negative value after offset: got {n}"
                        ),
                    ));
                }

                let q = n / m;
                let r = n % m;
                // Write q in unary (q zeros followed by a 1)
                for _ in 0..q {
                    core_data_writer.write_u32(0, 1)?;
                }
                core_data_writer.write_u32(1, 1)?;
                let b = 32 - ((*m - 1).leading_zeros());
                if b > 0 {
                    let threshold = (1i32 << b) - m;
                    if r < threshold {
                        core_data_writer.write_u32(r as u32, (b - 1) as usize)?;
                    } else {
                        let r2 = r + threshold;
                        core_data_writer.write_u32(r2 as u32, b as usize)?;
                    }
                }
                Ok(())
            }
            Self::GolombRice { offset, log2_m } => {
                if *log2_m < 0 || *log2_m >= 32 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid Golomb-Rice parameter: log2_m={log2_m}"),
                    ));
                }

                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let n = value + offset;

                if n < 0 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "GolombRice encoding requires non-negative value after offset: got {n}"
                        ),
                    ));
                }

                let log2_m = *log2_m as u32;
                let q = n >> log2_m;
                let r = n & ((1 << log2_m) - 1);
                // Write q in unary
                for _ in 0..q {
                    core_data_writer.write_u32(0, 1)?;
                }
                core_data_writer.write_u32(1, 1)?;
                core_data_writer.write_u32(r as u32, log2_m as usize)
            }
            Self::Subexp { offset, k } => {
                let value = i32::try_from(value)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                let n = value + offset;
                let k = *k as u32;

                if k >= 32 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("subexp encoding overflow: k={k}"),
                    ));
                }

                // Find the correct n_bits group
                let mut group = 0u32;
                let mut threshold = 1i32 << k;
                let mut prev_threshold = 0i32;
                loop {
                    if group < 2 {
                        if n < threshold {
                            // Write group 1-bits
                            for _ in 0..group {
                                core_data_writer.write_u32(1, 1)?;
                            }
                            core_data_writer.write_u32(0, 1)?;
                            let bits = k + group;
                            return core_data_writer
                                .write_u32((n - prev_threshold) as u32, bits as usize);
                        }
                        prev_threshold = threshold;
                        if group == 0 {
                            threshold = 1 << (k + 1);
                        } else {
                            threshold = prev_threshold + (1 << (k + 1));
                        }
                        group += 1;
                    } else {
                        let bits = k + group + group - 1;

                        if bits >= 32 {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("subexp encoding overflow: bits={bits}"),
                            ));
                        }

                        let max_val = prev_threshold + (1 << bits);
                        if n < max_val {
                            for _ in 0..group {
                                core_data_writer.write_u32(1, 1)?;
                            }
                            core_data_writer.write_u32(0, 1)?;
                            let extra = n - prev_threshold;
                            return core_data_writer.write_u32(extra as u32, bits as usize);
                        }
                        prev_threshold = max_val;
                        group += 1;
                    }
                }
            }
            Self::VarintUnsigned {
                block_content_id,
                offset,
            } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let adjusted = value - offset;
                let n = u64::try_from(adjusted)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                write_uint7_64(dst, n)
            }
            Self::VarintSigned {
                block_content_id,
                offset,
            } => {
                let dst = external_data_writers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                let adjusted = value - offset;
                write_sint7_64(dst, adjusted)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::container::compression_header::Encoding;

    #[test]
    fn test_decode() -> io::Result<()> {
        fn t(
            core_data: Option<&[u8]>,
            encoding: &Encoding<Integer>,
            expected: i64,
        ) -> io::Result<()> {
            let core_data = core_data.unwrap_or(&[0b10000000]);
            let mut core_data_reader = BitReader::new(core_data);

            let external_data = [0x0d];
            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(1, &external_data[..]);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        t(
            None,
            &Encoding::new(Integer::External {
                block_content_id: 1,
            }),
            13,
        )?;
        t(
            None,
            &Encoding::new(Integer::huffman(vec![0x4e], vec![0])),
            0x4e,
        )?;
        t(None, &Encoding::new(Integer::Beta { offset: 1, len: 3 }), 3)?;
        t(
            Some(&[0b00011010]),
            &Encoding::new(Integer::Gamma { offset: 5 }),
            8,
        )?;
        t(None, &Encoding::new(Integer::ConstInt { value: 42 }), 42)?;

        Ok(())
    }

    #[test]
    fn test_decode_varint_unsigned() -> io::Result<()> {
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        // uint7: 0x81 0x00 = 128
        let external_data = [0x81, 0x00];
        let mut external_data_readers = ExternalDataReaders::new();
        external_data_readers.insert(1, &external_data[..]);

        let encoding = Encoding::new(Integer::VarintUnsigned {
            block_content_id: 1,
            offset: 0,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 128);

        Ok(())
    }

    #[test]
    fn test_decode_varint_unsigned_with_offset() -> io::Result<()> {
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        // uint7: 0x05 = 5, offset = -1, result = 5 + (-1) = 4
        let external_data = [0x05];
        let mut external_data_readers = ExternalDataReaders::new();
        external_data_readers.insert(1, &external_data[..]);

        let encoding = Encoding::new(Integer::VarintUnsigned {
            block_content_id: 1,
            offset: -1,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 4);

        Ok(())
    }

    #[test]
    fn test_decode_varint_signed() -> io::Result<()> {
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        // sint7 zigzag: 0x03 = -2
        let external_data = [0x03];
        let mut external_data_readers = ExternalDataReaders::new();
        external_data_readers.insert(1, &external_data[..]);

        let encoding = Encoding::new(Integer::VarintSigned {
            block_content_id: 1,
            offset: 0,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, -2);

        Ok(())
    }

    #[test]
    fn test_decode_varint_signed_with_offset() -> io::Result<()> {
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        // sint7 zigzag: 0x04 = 2, offset = 10, result = 2 + 10 = 12
        let external_data = [0x04];
        let mut external_data_readers = ExternalDataReaders::new();
        external_data_readers.insert(1, &external_data[..]);

        let encoding = Encoding::new(Integer::VarintSigned {
            block_content_id: 1,
            offset: 10,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 12);

        Ok(())
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        fn t(
            encoding: &Encoding<Integer>,
            value: i64,
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
            &Encoding::new(Integer::External {
                block_content_id: 1,
            }),
            0x0d,
            &[],
            &[0x0d],
        )?;

        Ok(())
    }

    #[test]
    fn test_decode_golomb() -> io::Result<()> {
        // Golomb with m=5, offset=0
        // q=2 (two 0-bits then 1-bit), b=ceil(log2(5))=3
        // r from 2 bits: threshold = 2^3 - 5 = 3
        // bits: 00 1 01 => q=2, r=1 (1 < 3, so value = 2*5 + 1 = 11)
        let core_data = [0b00101000];
        let mut core_data_reader = BitReader::new(&core_data[..]);
        let mut external_data_readers = ExternalDataReaders::new();

        let encoding = Encoding::new(Integer::Golomb { offset: 0, m: 5 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 11);

        // Golomb with m=5, offset=10
        // Same bit pattern => value = 11 - 10 = 1
        let core_data = [0b00101000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Golomb { offset: 10, m: 5 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 1);

        // Golomb with m=5, offset=0, r >= threshold
        // bits: 1 11 0 => q=0, r from 2 bits = 3, 3 >= threshold(3),
        //   r = (3<<1)|0 = 6, value = 0*5 + 6 - 3 = 3
        let core_data = [0b11100000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Golomb { offset: 0, m: 5 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 3);

        Ok(())
    }

    #[test]
    fn test_decode_golomb_rice() -> io::Result<()> {
        // GolombRice with log2_m=3, offset=0
        // q from unary (count 0-bits): bits 0,0,1 => q=2
        // r from 3 bits: 101 = 5
        // value = (2 << 3) | 5 = 21
        let core_data = [0b00110100];
        let mut core_data_reader = BitReader::new(&core_data[..]);
        let mut external_data_readers = ExternalDataReaders::new();

        let encoding = Encoding::new(Integer::GolombRice {
            offset: 0,
            log2_m: 3,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 21);

        // GolombRice with log2_m=2, offset=5
        // bits: 1 10 => q=0, r=2, value = (0 << 2) | 2 - 5 = -3
        let core_data = [0b11000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::GolombRice {
            offset: 5,
            log2_m: 2,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, -3);

        Ok(())
    }

    #[test]
    fn test_decode_subexp() -> io::Result<()> {
        // Subexp with k=2, offset=0
        // n=0 (first bit is 0): read k+0=2 bits
        // bits: 0 11 => n=0, value = 3
        let core_data = [0b01100000];
        let mut core_data_reader = BitReader::new(&core_data[..]);
        let mut external_data_readers = ExternalDataReaders::new();

        let encoding = Encoding::new(Integer::Subexp { offset: 0, k: 2 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 3);

        // Subexp with k=2, offset=0, n=1
        // bits: 1 0 001 => n=1, read k+1=3 bits = 1, value = 1
        let core_data = [0b10001000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Subexp { offset: 0, k: 2 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 1);

        // Subexp with k=2, offset=3
        // bits: 0 10 => n=0, value from 2 bits = 2, result = 2 - 3 = -1
        let core_data = [0b01000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Subexp { offset: 3, k: 2 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, -1);

        // Subexp with k=0, offset=0
        // n=0 (first bit 0): read k+0=0 bits => value=0
        let core_data = [0b00000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Subexp { offset: 0, k: 0 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 0);

        // Subexp with k=0, offset=0, n=1
        // bits: 1 0 0 => n=1, read k+1=1 bit = 0, value = 0
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Subexp { offset: 0, k: 0 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 0);

        Ok(())
    }

    #[test]
    fn test_decode_golomb_with_m_1() -> io::Result<()> {
        // Golomb with m=1, offset=0 (degenerate case: unary coding)
        // With m=1: b=ceil(log2(1))=0, threshold=2^0-1=0
        // q from unary: count 0-bits before 1-bit
        // bits: 1 => q=0, no remainder bits, value = 0*1 + 0 = 0
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);
        let mut external_data_readers = ExternalDataReaders::new();

        let encoding = Encoding::new(Integer::Golomb { offset: 0, m: 1 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 0);

        // bits: 0 1 => q=1, value = 1*1 + 0 = 1
        let core_data = [0b01000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Golomb { offset: 0, m: 1 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 1);

        // bits: 0 0 0 1 => q=3, value = 3
        let core_data = [0b00010000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::Golomb { offset: 0, m: 1 });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 3);

        Ok(())
    }

    #[test]
    fn test_decode_golomb_rice_with_log2_m_0() -> io::Result<()> {
        // GolombRice with log2_m=0, offset=0 (equivalent to unary coding)
        // q from unary, r from 0 bits = 0
        // bits: 1 => q=0, value = (0 << 0) | 0 = 0
        let core_data = [0b10000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);
        let mut external_data_readers = ExternalDataReaders::new();

        let encoding = Encoding::new(Integer::GolombRice {
            offset: 0,
            log2_m: 0,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 0);

        // bits: 0 1 => q=1, value = (1 << 0) | 0 = 1
        let core_data = [0b01000000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::GolombRice {
            offset: 0,
            log2_m: 0,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, 1);

        // GolombRice with log2_m=0, offset=5
        // bits: 0 0 1 => q=2, value = 2 - 5 = -3
        let core_data = [0b00100000];
        let mut core_data_reader = BitReader::new(&core_data[..]);

        let encoding = Encoding::new(Integer::GolombRice {
            offset: 5,
            log2_m: 0,
        });
        let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;
        assert_eq!(actual, -3);

        Ok(())
    }
}
