use std::{io, num::NonZero};

use crate::codecs::rans_nx16::ALPHABET_SIZE;

// § 3.5 "rANS Nx16 Bit Packing" (2023-03-15): "It is not permitted to have `nsym` > 16 [...] as
// bit packing is not possible."
const MAX_SYMBOL_COUNT: usize = 16;

pub struct Context {
    pub symbol_count: NonZero<usize>,
    pub alphabet: [bool; ALPHABET_SIZE],
    pub mapping_table: [u8; ALPHABET_SIZE],
}

pub fn build_context(src: &[u8]) -> io::Result<Context> {
    let alphabet = build_alphabet(src);
    let (mapping_table, symbol_count) = build_mapping_table(&alphabet)?;

    Ok(Context {
        symbol_count,
        alphabet,
        mapping_table,
    })
}

fn build_alphabet(src: &[u8]) -> [bool; ALPHABET_SIZE] {
    let mut alphabet = [false; ALPHABET_SIZE];

    for &sym in src {
        alphabet[usize::from(sym)] = true;
    }

    alphabet
}

fn build_mapping_table(
    alphabet: &[bool; ALPHABET_SIZE],
) -> io::Result<([u8; ALPHABET_SIZE], NonZero<usize>)> {
    let mut mapping_table = [0; ALPHABET_SIZE];
    let mut symbol_count = 0;

    for (sym, _) in alphabet.iter().enumerate().filter(|(_, a)| **a) {
        mapping_table[sym] = symbol_count;
        symbol_count += 1;

        let n = usize::from(symbol_count);

        if n > MAX_SYMBOL_COUNT {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }
    }

    let n = NonZero::new(usize::from(symbol_count))
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidInput))?;

    Ok((mapping_table, n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_context() -> io::Result<()> {
        let src = b"ndls";

        let mut alphabet = [false; ALPHABET_SIZE];
        alphabet[usize::from(b'n')] = true;
        alphabet[usize::from(b'd')] = true;
        alphabet[usize::from(b'l')] = true;
        alphabet[usize::from(b's')] = true;

        let mut mapping_table = [0; ALPHABET_SIZE];
        mapping_table[usize::from(b'd')] = 0;
        mapping_table[usize::from(b'l')] = 1;
        mapping_table[usize::from(b'n')] = 2;
        mapping_table[usize::from(b's')] = 3;

        let expected = Context {
            symbol_count: const { NonZero::new(4).unwrap() },
            alphabet,
            mapping_table,
        };

        let actual = build_context(src)?;
        assert_eq!(actual.symbol_count, expected.symbol_count);
        assert_eq!(actual.alphabet, expected.alphabet);
        assert_eq!(actual.mapping_table, expected.mapping_table);

        Ok(())
    }

    #[test]
    fn test_build_alphabet() {
        let src = [];
        let expected = [false; ALPHABET_SIZE];
        assert_eq!(build_alphabet(&src), expected);

        let src = b"ndls";
        let mut expected = [false; ALPHABET_SIZE];
        expected[usize::from(b'n')] = true;
        expected[usize::from(b'd')] = true;
        expected[usize::from(b'l')] = true;
        expected[usize::from(b's')] = true;
        assert_eq!(build_alphabet(src), expected);
    }

    #[test]
    fn test_build_mapping_table() -> io::Result<()> {
        let src = b"ndls";
        let alphabet = build_alphabet(src);

        let mut expected_mapping_table = [0; ALPHABET_SIZE];
        expected_mapping_table[usize::from(b'd')] = 0;
        expected_mapping_table[usize::from(b'l')] = 1;
        expected_mapping_table[usize::from(b'n')] = 2;
        expected_mapping_table[usize::from(b's')] = 3;

        let expected_symbol_count = const { NonZero::new(4).unwrap() };

        let expected = (expected_mapping_table, expected_symbol_count);
        assert_eq!(build_mapping_table(&alphabet)?, expected);

        // 0 symbols
        let src = [];
        let alphabet = build_alphabet(&src);
        assert!(matches!(
            build_mapping_table(&alphabet),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        // > 16 symbols
        let src = b"abcdefghijklmnopq";
        let alphabet = build_alphabet(src);
        assert!(matches!(
            build_mapping_table(&alphabet),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
