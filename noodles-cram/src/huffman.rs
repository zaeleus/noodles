use std::{
    collections::HashMap,
    io::{self, Read},
};

use crate::{num::Itf8, BitReader};

type CodeBook = HashMap<Itf8, (Itf8, usize)>;

pub struct CanonicalHuffmanDecoder {
    code_book: CodeBook,
}

impl CanonicalHuffmanDecoder {
    pub fn new(alphabet: &[Itf8], bit_lens: &[Itf8]) -> Self {
        let code_book = build_canonical_code_book(alphabet, bit_lens);
        Self { code_book }
    }

    pub fn read<R>(&self, reader: &mut BitReader<R>) -> io::Result<Itf8>
    where
        R: Read,
    {
        let mut code_book_by_len = HashMap::new();

        for (symbol, (code, len)) in &self.code_book {
            let entries = code_book_by_len.entry(len).or_insert(Vec::new());
            entries.push((symbol, code));
        }

        let sorted_lens: Vec<_> = code_book_by_len.keys().cloned().collect();

        let mut prev_len = 0;
        let mut input_code = 0;

        for &len in sorted_lens {
            input_code <<= len - prev_len;

            let b = reader.read_u32(len)? as i32;
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

fn build_canonical_code_book(alphabet: &[Itf8], bit_lens: &[Itf8]) -> CodeBook {
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

        code_book.insert(symbol, (code, bit_len as usize));

        code += 1;
        prev_bit_len = bit_len;
    }

    code_book
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_canonical_code_book() {
        let symbols = [65, 66, 67, 68, 69, 70];
        let bit_lens = [1, 3, 3, 3, 4, 4];

        let code_book = build_canonical_code_book(&symbols, &bit_lens);

        assert_eq!(code_book.len(), 6);

        // assert_eq!(code_book[&65], 0b0);
        // assert_eq!(code_book[&66], 0b100);
        // assert_eq!(code_book[&67], 0b101);
        // assert_eq!(code_book[&68], 0b110);
        // assert_eq!(code_book[&69], 0b1110);
        // assert_eq!(code_book[&70], 0b1111);
    }
}
