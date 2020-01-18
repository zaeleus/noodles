use std::ops::{Deref, Index};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    Eq,
    A,
    C,
    M,
    G,
    R,
    S,
    V,
    T,
    W,
    Y,
    H,
    K,
    D,
    B,
    N,
}

impl Base {
    // https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
    pub fn complement(self) -> Self {
        match self {
            Base::Eq => Base::Eq,
            Base::A => Base::T,
            Base::C => Base::G,
            Base::G => Base::C,
            Base::T => Base::A,
            Base::W => Base::W,
            Base::S => Base::S,
            Base::M => Base::K,
            Base::K => Base::M,
            Base::R => Base::Y,
            Base::Y => Base::R,
            Base::B => Base::V,
            Base::D => Base::H,
            Base::H => Base::D,
            Base::V => Base::B,
            Base::N => Base::N,
        }
    }

    pub fn symbol(self) -> char {
        match self {
            Base::Eq => '=',
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::W => 'W',
            Base::S => 'S',
            Base::M => 'M',
            Base::K => 'K',
            Base::R => 'R',
            Base::Y => 'Y',
            Base::B => 'B',
            Base::D => 'D',
            Base::H => 'H',
            Base::V => 'V',
            Base::N => 'N',
        }
    }
}

impl Into<char> for Base {
    fn into(self) -> char {
        self.symbol()
    }
}

static BASES: &[Base] = &[
    Base::Eq,
    Base::A,
    Base::C,
    Base::M,
    Base::G,
    Base::R,
    Base::S,
    Base::V,
    Base::T,
    Base::W,
    Base::Y,
    Base::H,
    Base::K,
    Base::D,
    Base::B,
    Base::N,
];

#[derive(Debug)]
pub struct Sequence<'a> {
    seq: &'a [u8],
    n_chars: usize,
}

impl<'a> Sequence<'a> {
    pub fn new(seq: &[u8], n_chars: usize) -> Sequence {
        Sequence { seq, n_chars }
    }

    pub fn n_chars(&self) -> usize {
        self.n_chars
    }

    pub fn bases(&self) -> Bases {
        Bases {
            sequence: self,
            head: 0,
            tail: self.n_chars - 1,
            remaining: self.n_chars,
        }
    }

    pub fn symbols(&self) -> Symbols<Bases> {
        Symbols { iter: self.bases() }
    }
}

impl<'a> Deref for Sequence<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.seq
    }
}

pub struct Bases<'a> {
    sequence: &'a Sequence<'a>,
    head: usize,
    tail: usize,
    remaining: usize,
}

impl<'a> Iterator for Bases<'a> {
    type Item = Base;

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.head];
        self.head += 1;
        self.remaining -= 1;
        Some(symbol)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}

impl<'a> DoubleEndedIterator for Bases<'a> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        let symbol = self.sequence[self.tail];
        self.tail -= 1;
        self.remaining -= 1;
        Some(symbol)
    }
}

pub struct Symbols<I> {
    iter: I,
}

impl<I: Iterator<Item = Base>> Iterator for Symbols<I> {
    type Item = char;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|b| b.symbol())
    }
}

pub struct Complement<I> {
    iter: I,
}

impl<I: Iterator<Item = Base>> Complement<I> {
    pub fn new(iter: I) -> Complement<I> {
        Complement { iter }
    }
}

impl<I: Iterator<Item = Base>> Iterator for Complement<I> {
    type Item = Base;

    fn next(&mut self) -> Option<Base> {
        self.iter.next().map(|b| b.complement())
    }
}

impl<'a> Index<usize> for Sequence<'a> {
    type Output = Base;

    fn index(&self, i: usize) -> &Self::Output {
        let j = i / 2;
        let b = self.seq[j];

        let k = if i % 2 == 0 {
            (b & 0xf0) >> 4
        } else {
            b & 0x0f
        };

        &BASES[k as usize]
    }
}

#[cfg(test)]
mod base_tests {
    use super::Base;

    #[test]
    fn test_complement() {
        assert_eq!(Base::Eq.complement(), Base::Eq);
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::W.complement(), Base::W);
        assert_eq!(Base::S.complement(), Base::S);
        assert_eq!(Base::M.complement(), Base::K);
        assert_eq!(Base::K.complement(), Base::M);
        assert_eq!(Base::R.complement(), Base::Y);
        assert_eq!(Base::Y.complement(), Base::R);
        assert_eq!(Base::B.complement(), Base::V);
        assert_eq!(Base::D.complement(), Base::H);
        assert_eq!(Base::H.complement(), Base::D);
        assert_eq!(Base::V.complement(), Base::B);
        assert_eq!(Base::N.complement(), Base::N);
    }

    #[test]
    fn test_symbol() {
        assert_eq!(Base::Eq.symbol(), '=');
        assert_eq!(Base::A.symbol(), 'A');
        assert_eq!(Base::C.symbol(), 'C');
        assert_eq!(Base::G.symbol(), 'G');
        assert_eq!(Base::T.symbol(), 'T');
        assert_eq!(Base::W.symbol(), 'W');
        assert_eq!(Base::S.symbol(), 'S');
        assert_eq!(Base::M.symbol(), 'M');
        assert_eq!(Base::K.symbol(), 'K');
        assert_eq!(Base::R.symbol(), 'R');
        assert_eq!(Base::Y.symbol(), 'Y');
        assert_eq!(Base::B.symbol(), 'B');
        assert_eq!(Base::D.symbol(), 'D');
        assert_eq!(Base::H.symbol(), 'H');
        assert_eq!(Base::V.symbol(), 'V');
        assert_eq!(Base::N.symbol(), 'N');
    }
}

#[cfg(test)]
mod sequence_tests {
    use super::*;

    #[test]
    fn test_bases() {
        let data = [0x18, 0x42];
        let sequence = Sequence::new(&data, 4);

        let mut bases = sequence.bases();

        assert_eq!(bases.next(), Some(Base::A));
        assert_eq!(bases.next(), Some(Base::T));
        assert_eq!(bases.next(), Some(Base::G));
        assert_eq!(bases.next(), Some(Base::C));
        assert_eq!(bases.next(), None);
    }

    #[test]
    fn test_symbols() {
        let data = [0x18, 0x42];
        let sequence = Sequence::new(&data, 4);

        let mut bases = sequence.symbols();

        assert_eq!(bases.next(), Some('A'));
        assert_eq!(bases.next(), Some('T'));
        assert_eq!(bases.next(), Some('G'));
        assert_eq!(bases.next(), Some('C'));
        assert_eq!(bases.next(), None);
    }
}
