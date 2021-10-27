//! BAM record sequence and bases.

mod base;
mod bases;

pub use self::{base::Base, bases::Bases};

use std::{
    fmt,
    ops::{Deref, DerefMut},
};

use noodles_sam as sam;

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

/// A BAM record sequence.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Sequence {
    seq: Vec<u8>,
    base_count: usize,
}

impl Sequence {
    /// Creates a sequence by wrapping raw sequence data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    ///
    /// let data = vec![0x12, 0x48]; // ACGT
    /// let sequence = Sequence::new(data, 4);
    /// assert_eq!(**sequence, [0x12, 0x48]);
    /// ```
    pub fn new(seq: Vec<u8>, base_count: usize) -> Self {
        Self { seq, base_count }
    }

    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    ///
    /// let data = vec![0x12, 0x48]; // ACGT
    /// let sequence = Sequence::new(data, 4);
    ///
    /// assert_eq!(sequence.base_count(), 4);
    /// ```
    pub fn base_count(&self) -> usize {
        self.base_count
    }

    /// Returns a mutable reference to the base count.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    /// let mut sequence = Sequence::default();
    /// sequence.set_base_count(2);
    /// assert_eq!(sequence.base_count(), 2);
    /// ```
    pub fn set_base_count(&mut self, base_count: usize) {
        self.base_count = base_count;
    }

    /// Returns a reference to the base at the given index.
    ///
    /// If the index is out of bounds, this returns `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::{sequence::Base, Sequence};
    ///
    /// let data = vec![0x12, 0x48]; // ACGT
    /// let sequence = Sequence::new(data, 4);
    ///
    /// assert_eq!(sequence.get(1), Some(&Base::C));
    /// assert_eq!(sequence.get(8), None);
    /// ```
    pub fn get(&self, i: usize) -> Option<&Base> {
        let j = i / 2;
        let b = self.seq.get(j)?;

        if i % 2 == 0 {
            let k = (b & 0xf0) >> 4;
            Some(&BASES[k as usize])
        } else if i < self.base_count {
            let k = b & 0x0f;
            Some(&BASES[k as usize])
        } else {
            None
        }
    }

    /// Returns a iterator over the bases in this sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::{sequence::Base, Sequence};
    ///
    /// let data = vec![0x12, 0x48]; // ACGT
    /// let sequence = Sequence::new(data, 4);
    ///
    /// let mut bases = sequence.bases();
    ///
    /// assert_eq!(bases.next(), Some(Base::A));
    /// assert_eq!(bases.next(), Some(Base::C));
    /// assert_eq!(bases.next(), Some(Base::G));
    /// assert_eq!(bases.next(), Some(Base::T));
    /// assert_eq!(bases.next(), None);
    /// ```
    pub fn bases(&self) -> Bases<'_> {
        Bases::new(self)
    }
}

impl<'a> fmt::Debug for Sequence {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_list().entries(self.bases()).finish()
    }
}

impl Deref for Sequence {
    type Target = Vec<u8>;

    fn deref(&self) -> &Self::Target {
        &self.seq
    }
}

impl DerefMut for Sequence {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.seq
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in self.bases() {
            write!(f, "{}", base)?;
        }

        Ok(())
    }
}

impl From<&Sequence> for sam::record::Sequence {
    fn from(sequence: &Sequence) -> Self {
        let sam_bases: Vec<_> = sequence.bases().map(|b| b.into()).collect();
        Self::from(sam_bases)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get() {
        let data = vec![0x18, 0x40];
        let sequence = Sequence::new(data, 3);
        assert_eq!(sequence.get(0), Some(&Base::A));
        assert_eq!(sequence.get(1), Some(&Base::T));
        assert_eq!(sequence.get(2), Some(&Base::G));
        assert_eq!(sequence.get(3), None);
        assert_eq!(sequence.get(8), None);
    }

    #[test]
    fn test_bases() {
        let data = vec![0x18, 0x42];
        let sequence = Sequence::new(data, 4);

        let mut bases = sequence.bases();

        assert_eq!(bases.next(), Some(Base::A));
        assert_eq!(bases.next(), Some(Base::T));
        assert_eq!(bases.next(), Some(Base::G));
        assert_eq!(bases.next(), Some(Base::C));
        assert_eq!(bases.next(), None);
    }

    #[test]
    fn test_bases_with_empty_sequence() {
        let sequence = Sequence::new(Vec::new(), 0);
        assert!(sequence.bases().next().is_none());
    }

    #[test]
    fn test_fmt() {
        let data = vec![0x18, 0x42];
        let sequence = Sequence::new(data, 4);
        assert_eq!(sequence.to_string(), "ATGC");
    }

    #[test]
    fn test_from_sequence_for_sam_record_sequence() {
        use sam::record::{sequence::Base as SamBase, Sequence as SamSequence};

        // ATGC
        let data = vec![0x18, 0x42];
        let sequence = Sequence::new(data, 4);

        let actual = SamSequence::from(&sequence);
        let expected = SamSequence::from(vec![SamBase::A, SamBase::T, SamBase::G, SamBase::C]);

        assert_eq!(actual, expected);
    }
}
