use std::{
    io,
    num::{NonZero, Saturating},
};

use crate::{codecs::rans_nx16::ALPHABET_SIZE, io::writer::num::write_u8};

pub struct Context {
    pub alphabet: [i32; ALPHABET_SIZE],
    pub dst: Vec<u8>,
}

pub fn build_context(src: &[u8]) -> io::Result<Context> {
    let mut alphabet = [Saturating::default(); ALPHABET_SIZE];

    for window in src.windows(2) {
        let (prev_sym, sym) = (window[0], window[1]);
        let delta = if sym == prev_sym { 1 } else { -1 };
        alphabet[usize::from(sym)] += delta;
    }

    let alphabet = alphabet.map(|n| n.0);

    let n = alphabet.iter().filter(|&&n| n > 0).count();
    let symbol_count =
        NonZero::new(n).ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

    let mut dst = Vec::new();
    write_alphabet(&mut dst, &alphabet, symbol_count)?;

    Ok(Context { alphabet, dst })
}

fn write_alphabet(
    dst: &mut Vec<u8>,
    alphabet: &[i32; ALPHABET_SIZE],
    symbol_count: NonZero<usize>,
) -> io::Result<()> {
    write_symbol_count(dst, symbol_count)?;

    for (sym, _) in alphabet.iter().enumerate().filter(|(_, n)| **n > 0) {
        // SAFETY: `sym` < `ALPHABET_SIZE`.
        write_u8(dst, sym as u8)?;
    }

    Ok(())
}

fn write_symbol_count(dst: &mut Vec<u8>, symbol_count: NonZero<usize>) -> io::Result<()> {
    let mut n = symbol_count.get();

    if n >= ALPHABET_SIZE {
        n = 0;
    }

    // SAFETY: `n` < `ALPHABET_SIZE`.
    write_u8(dst, n as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_context() -> io::Result<()> {
        let src = b"nnndlllllss";
        let ctx = build_context(src)?;

        let mut expected_alphabet = [0; ALPHABET_SIZE];
        expected_alphabet[usize::from(b'd')] = -1;
        expected_alphabet[usize::from(b'l')] = 3;
        expected_alphabet[usize::from(b'n')] = 2;
        expected_alphabet[usize::from(b's')] = 0;
        assert_eq!(ctx.alphabet, expected_alphabet);

        assert_eq!(
            ctx.dst,
            [
                0x02, // symbol count
                b'l', b'n', // alphabet
            ]
        );

        Ok(())
    }

    #[test]
    fn test_write_alphabet() -> io::Result<()> {
        let mut alphabet = [0; ALPHABET_SIZE];
        alphabet[usize::from(b'd')] = -1;
        alphabet[usize::from(b'n')] = 1;

        let symbol_count = const { NonZero::new(1).unwrap() };

        let mut dst = Vec::new();
        write_alphabet(&mut dst, &alphabet, symbol_count)?;

        assert_eq!(
            dst,
            [
                0x01, // symbol count = 1
                b'n', // alphabet
            ]
        );

        Ok(())
    }

    #[test]
    fn test_write_symbol_count() -> io::Result<()> {
        let mut dst = Vec::new();

        dst.clear();
        write_symbol_count(&mut dst, const { NonZero::new(1).unwrap() })?;
        assert_eq!(dst, [0x01]);

        dst.clear();
        write_symbol_count(&mut dst, const { NonZero::new(256).unwrap() })?;
        assert_eq!(dst, [0x00]);

        Ok(())
    }
}
