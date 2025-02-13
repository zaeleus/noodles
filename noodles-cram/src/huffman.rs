use std::{collections::HashMap, io};

use crate::io::BitReader;

type CodeBook = HashMap<i32, (i32, u32)>;

pub struct CanonicalHuffmanDecoder {
    code_book: CodeBook,
}

impl CanonicalHuffmanDecoder {
    pub fn new(alphabet: &[i32], bit_lens: &[u32]) -> Self {
        let code_book = build_canonical_code_book(alphabet, bit_lens);
        Self { code_book }
    }

    pub fn decode(&self, reader: &mut BitReader<'_>) -> io::Result<i32> {
        let mut code_book_by_len = HashMap::new();

        for (symbol, (code, len)) in &self.code_book {
            let entries = code_book_by_len.entry(len).or_insert_with(Vec::new);
            entries.push((symbol, code));
        }

        let sorted_lens = {
            let mut lens: Vec<_> = code_book_by_len.keys().copied().collect();
            lens.sort_unstable();
            lens
        };

        let mut prev_len = 0;
        let mut input_code = 0;

        for &len in sorted_lens {
            input_code <<= len - prev_len;

            let b = reader.read_u32(len - prev_len)? as i32;
            input_code |= b;

            let entry = code_book_by_len[&len]
                .iter()
                .find(|(_, code)| input_code == **code);

            if let Some((symbol, _)) = entry {
                return Ok(**symbol);
            }

            prev_len = len;
        }

        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "could not find symbol",
        ))
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
