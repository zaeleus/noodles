use std::{collections::HashMap, io};

use crate::io::{BitReader, BitWriter};

type CodeBook = HashMap<i32, (i32, u32)>;

#[derive(Clone, Debug)]
pub struct CanonicalHuffmanDecoder {
    // Pre-computed decode table: vec of (bit_len, vec of (code, symbol)), sorted by bit_len.
    decode_table: Vec<(u32, Vec<(i32, i32)>)>,
}

impl CanonicalHuffmanDecoder {
    pub fn new(alphabet: &[i32], bit_lens: &[u32]) -> Self {
        let code_book = build_canonical_code_book(alphabet, bit_lens);

        let mut by_len: HashMap<u32, Vec<(i32, i32)>> = HashMap::new();

        for (symbol, (code, len)) in &code_book {
            by_len.entry(*len).or_default().push((*code, *symbol));
        }

        let mut decode_table: Vec<_> = by_len.into_iter().collect();
        decode_table.sort_unstable_by_key(|(len, _)| *len);

        Self { decode_table }
    }

    pub fn decode(&self, reader: &mut BitReader<'_>) -> io::Result<i32> {
        let mut prev_len = 0;
        let mut input_code = 0;

        for (len, entries) in &self.decode_table {
            input_code <<= len - prev_len;

            let b = reader.read_i32(len - prev_len)?;
            input_code |= b;

            if let Some((_, symbol)) = entries.iter().find(|(code, _)| input_code == *code) {
                return Ok(*symbol);
            }

            prev_len = *len;
        }

        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "could not find symbol",
        ))
    }
}

#[derive(Clone, Debug)]
pub struct CanonicalHuffmanEncoder {
    code_book: CodeBook,
}

impl CanonicalHuffmanEncoder {
    pub fn new(alphabet: &[i32], bit_lens: &[u32]) -> Self {
        let code_book = build_canonical_code_book(alphabet, bit_lens);
        Self { code_book }
    }

    pub fn encode(&self, writer: &mut BitWriter, value: i32) -> io::Result<()> {
        let &(code, bit_len) = self.code_book.get(&value).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("symbol not in code book: {value}"),
            )
        })?;

        writer.write_u32(code as u32, bit_len as usize)
    }
}

fn build_canonical_code_book(alphabet: &[i32], bit_lens: &[u32]) -> CodeBook {
    let sorted_alphabet = {
        let mut pairs: Vec<_> = alphabet.iter().zip(bit_lens.iter()).collect();
        pairs.sort_by_key(|&(symbol, bit_len)| (bit_len, symbol));
        pairs
    };

    let mut code_book = CodeBook::with_capacity(sorted_alphabet.len());

    let mut code = 0;
    let mut prev_bit_len = *sorted_alphabet[0].1;

    for (&symbol, &bit_len) in sorted_alphabet {
        if bit_len > prev_bit_len {
            code <<= bit_len - prev_bit_len;
        }

        code_book.insert(symbol, (code, bit_len));

        code += 1;
        prev_bit_len = bit_len;
    }

    code_book
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        let symbols = [0x4e, 0x44, 0x4c];
        let bit_lens = [1, 2, 2];
        let decoder = CanonicalHuffmanDecoder::new(&symbols, &bit_lens);

        let data = [0b01011000];
        let mut reader = BitReader::new(&data[..]);

        assert_eq!(decoder.decode(&mut reader)?, 0x4e);
        assert_eq!(decoder.decode(&mut reader)?, 0x44);
        assert_eq!(decoder.decode(&mut reader)?, 0x4c);
        assert_eq!(decoder.decode(&mut reader)?, 0x4e);

        Ok(())
    }

    #[test]
    fn test_encode() -> io::Result<()> {
        let symbols = [0x4e, 0x44, 0x4c];
        let bit_lens = [1, 2, 2];
        let encoder = CanonicalHuffmanEncoder::new(&symbols, &bit_lens);

        let mut writer = BitWriter::default();
        encoder.encode(&mut writer, 0x4e)?;
        encoder.encode(&mut writer, 0x44)?;
        encoder.encode(&mut writer, 0x4c)?;
        encoder.encode(&mut writer, 0x4e)?;

        let data = writer.finish()?;
        assert_eq!(data, [0b01011000]);

        Ok(())
    }

    #[test]
    fn test_round_trip() -> io::Result<()> {
        let symbols = [65, 66, 67, 68, 69, 70];
        let bit_lens = [1, 3, 3, 3, 4, 4];

        let encoder = CanonicalHuffmanEncoder::new(&symbols, &bit_lens);
        let decoder = CanonicalHuffmanDecoder::new(&symbols, &bit_lens);

        let values = [65, 66, 67, 68, 69, 70, 65, 65];

        let mut writer = BitWriter::default();
        for &v in &values {
            encoder.encode(&mut writer, v)?;
        }
        let data = writer.finish()?;

        let mut reader = BitReader::new(&data);
        for &expected in &values {
            let actual = decoder.decode(&mut reader)?;
            assert_eq!(actual, expected);
        }

        Ok(())
    }

    #[test]
    fn test_single_symbol_encode() -> io::Result<()> {
        let symbols = [42];
        let bit_lens = [0];
        let encoder = CanonicalHuffmanEncoder::new(&symbols, &bit_lens);

        let mut writer = BitWriter::default();
        encoder.encode(&mut writer, 42)?;
        encoder.encode(&mut writer, 42)?;

        let data = writer.finish()?;
        assert!(data.is_empty());

        Ok(())
    }

    #[test]
    fn test_build_canonical_code_book() {
        let symbols = [65, 66, 67, 68, 69, 70];
        let bit_lens = [1, 3, 3, 3, 4, 4];

        let code_book = build_canonical_code_book(&symbols, &bit_lens);

        assert_eq!(code_book.len(), 6);

        assert_eq!(code_book[&65], (0b0, 1));
        assert_eq!(code_book[&66], (0b100, 3));
        assert_eq!(code_book[&67], (0b101, 3));
        assert_eq!(code_book[&68], (0b110, 3));
        assert_eq!(code_book[&69], (0b1110, 4));
        assert_eq!(code_book[&70], (0b1111, 4));
    }
}
