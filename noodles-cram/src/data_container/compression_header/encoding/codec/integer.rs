use std::io;

use bytes::Buf;

use crate::{
    container::block,
    data_container::compression_header::encoding::Decode,
    huffman::CanonicalHuffmanDecoder,
    io::BitReader,
    reader::{num::get_itf8, record::ExternalDataReaders},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Integer {
    // block_content_id
    External(block::ContentId),
    // offset, m
    Golomb(i32, i32),
    // alphabet, bit_lens
    Huffman(Vec<i32>, Vec<u32>),
    // offset, len
    Beta(i32, u32),
    // offset, k
    Subexp(i32, i32),
    // offset, log2_m
    GolombRice(i32, i32),
    // offset
    Gamma(i32),
}

impl Decode for Integer {
    type Value = i32;

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
            Integer::External(block_content_id) => {
                let src = external_data_readers
                    .get_mut(block_content_id)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("missing external block: {block_content_id}"),
                        )
                    })?;

                get_itf8(src)
            }
            Integer::Huffman(alphabet, bit_lens) => {
                if alphabet.len() == 1 {
                    Ok(alphabet[0])
                } else {
                    let decoder = CanonicalHuffmanDecoder::new(alphabet, bit_lens);
                    decoder.decode(core_data_reader)
                }
            }
            Integer::Beta(offset, len) => {
                core_data_reader.read_u32(*len).map(|i| (i as i32 - offset))
            }
            Integer::Gamma(offset) => {
                let mut n = 0;

                while core_data_reader.read_bit()? == 0 {
                    n += 1;
                }

                let m = core_data_reader.read_u32(n)? as i32;
                let x = (1 << n) + m;

                Ok(x - offset)
            }
            _ => todo!("decode_itf8: {:?}", self),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_container::compression_header::Encoding;

    #[test]
    fn test_decode_itf8() -> io::Result<()> {
        fn t(
            core_data: Option<&[u8]>,
            encoding: &Encoding<Integer>,
            expected: i32,
        ) -> io::Result<()> {
            let core_data = core_data.unwrap_or(&[0b10000000]);
            let mut core_data_reader = BitReader::new(core_data);

            let external_data = [0x0d];
            let mut external_data_readers = ExternalDataReaders::new();
            external_data_readers.insert(block::ContentId::from(1), &external_data[..]);

            let actual = encoding.decode(&mut core_data_reader, &mut external_data_readers)?;

            assert_eq!(expected, actual);

            Ok(())
        }

        t(
            None,
            &Encoding::new(Integer::External(block::ContentId::from(1))),
            13,
        )?;
        t(
            None,
            &Encoding::new(Integer::Huffman(vec![0x4e], vec![0])),
            0x4e,
        )?;
        t(None, &Encoding::new(Integer::Beta(1, 3)), 3)?;
        t(Some(&[0b00011010]), &Encoding::new(Integer::Gamma(5)), 8)?;

        Ok(())
    }
}
