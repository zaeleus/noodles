mod base;
mod bases;

pub use self::{base::Base, bases::Bases};

use std::{
    fmt,
    ops::{Deref, Index},
};

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

pub struct Sequence<'a> {
    seq: &'a [u8],
    n_chars: usize,
}

impl<'a> Sequence<'a> {
    pub fn new(seq: &'a [u8], n_chars: usize) -> Self {
        Self { seq, n_chars }
    }

    pub fn n_chars(&self) -> usize {
        self.n_chars
    }

    pub fn get(&self, i: usize) -> Option<&Base> {
        let j = i / 2;
        let b = self.seq.get(j)?;

        let k = if i % 2 == 0 {
            (b & 0xf0) >> 4
        } else {
            b & 0x0f
        };

        BASES.get(k as usize)
    }

    pub fn bases(&self) -> Bases<'_> {
        Bases::new(self)
    }
}

impl<'a> fmt::Debug for Sequence<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_list().entries(self.bases()).finish()
    }
}

impl<'a> Deref for Sequence<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.seq
    }
}

impl<'a> fmt::Display for Sequence<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self.bases() {
            write!(f, "{}", base)?;
        }

        Ok(())
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
mod tests {
    use super::*;

    #[test]
    fn test_get() {
        let data = [0x18, 0x42];
        let sequence = Sequence::new(&data, 4);
        assert_eq!(sequence.get(2), Some(&Base::G));
        assert_eq!(sequence.get(8), None);
    }

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
    fn test_fmt() {
        let data = [0x18, 0x42];
        let sequence = Sequence::new(&data, 4);
        assert_eq!(sequence.to_string(), "ATGC");
    }
}
