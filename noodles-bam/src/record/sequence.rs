//! BAM record sequence and bases.

mod base;
mod bases;

pub use self::{base::Base, bases::Bases};

use std::{error, fmt};

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
    len: usize,
}

impl Sequence {
    /// Creates a sequence by wrapping raw sequence data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
    /// assert_eq!(sequence.as_ref(), [0x12, 0x48]);
    /// ```
    pub fn new(seq: Vec<u8>, len: usize) -> Self {
        Self { seq, len }
    }

    /// Returns whether the sequence contains any bases.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    ///
    /// let sequence = Sequence::default();
    /// assert!(sequence.is_empty());
    ///
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
    /// assert!(!sequence.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
    /// assert_eq!(sequence.len(), 4);
    /// ```
    #[deprecated(since = "0.8.0", note = "Use `Sequence::len` instead.")]
    pub fn base_count(&self) -> usize {
        self.len()
    }

    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
    /// assert_eq!(sequence.len(), 4);
    /// ```
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns a mutable reference to the base count.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Sequence;
    /// let mut sequence = Sequence::default();
    /// sequence.set_len(2);
    /// assert_eq!(sequence.len(), 2);
    /// ```
    pub fn set_len(&mut self, len: usize) {
        self.len = len;
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
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
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
        } else if i < self.len {
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
    /// let sequence = Sequence::new(vec![0x12, 0x48], 4); // ACGT
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

    /// Appends a base to the end of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::{sequence::Base, Sequence};
    ///
    /// let mut sequence = Sequence::new(vec![0x12, 0x40], 3); // ACG
    ///
    /// sequence.push(Base::T);
    /// assert_eq!(sequence.as_ref(), [0x12, 0x48]); // ACGT
    /// assert_eq!(sequence.len(), 4);
    ///
    /// sequence.push(Base::A);
    /// assert_eq!(sequence.as_ref(), [0x12, 0x48, 0x10]); //ACGTA
    /// assert_eq!(sequence.len(), 5);
    /// ```
    pub fn push(&mut self, base: Base) {
        let b = u8::from(base);

        if self.len % 2 == 0 {
            self.seq.push(b << 4);
        } else if let Some(l) = self.seq.last_mut() {
            *l |= b;
        }

        self.len += 1;
    }
}

impl<'a> fmt::Debug for Sequence {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_list().entries(self.bases()).finish()
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.seq
    }
}

impl AsMut<Vec<u8>> for Sequence {
    fn as_mut(&mut self) -> &mut Vec<u8> {
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

impl From<Vec<Base>> for Sequence {
    fn from(bases: Vec<Base>) -> Self {
        let data = bases
            .chunks(2)
            .map(|chunk| {
                let l = chunk[0];
                let r = chunk.get(1).copied().unwrap_or(Base::Eq);
                u8::from(l) << 4 | u8::from(r)
            })
            .collect();

        Self::new(data, bases.len())
    }
}

/// An error returned when a raw BAM record sequence fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// A base is invalid.
    InvalidBase(sam::record::sequence::base::TryFromCharError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBase(e) => write!(f, "invalid base: {}", e),
        }
    }
}

impl TryFrom<&[u8]> for Sequence {
    type Error = ParseError;

    fn try_from(data: &[u8]) -> Result<Self, Self::Error> {
        data.iter()
            .map(|&b| {
                sam::record::sequence::Base::try_from(char::from(b))
                    .map(Base::from)
                    .map_err(ParseError::InvalidBase)
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Self::from)
    }
}

impl From<&Sequence> for sam::record::Sequence {
    fn from(sequence: &Sequence) -> Self {
        let bases: Vec<_> = sequence.bases().map(|b| b.into()).collect();
        Self::from(bases)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get() {
        let sequence = Sequence::new(vec![0x18, 0x40], 3); // ATG

        assert_eq!(sequence.get(0), Some(&Base::A));
        assert_eq!(sequence.get(1), Some(&Base::T));
        assert_eq!(sequence.get(2), Some(&Base::G));
        assert_eq!(sequence.get(3), None);
        assert_eq!(sequence.get(8), None);
    }

    #[test]
    fn test_bases() {
        let sequence = Sequence::new(vec![0x18, 0x42], 4); // ATGC

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
        let sequence = Sequence::new(vec![0x18, 0x42], 4); // ATGC
        assert_eq!(sequence.to_string(), "ATGC");
    }

    #[test]
    fn test_from_vec_base_for_sequence() {
        let actual = Sequence::from(vec![Base::A, Base::T, Base::G]);
        let expected = Sequence::new(vec![0x18, 0x40], 3);
        assert_eq!(actual, expected);

        let actual = Sequence::from(vec![Base::A, Base::T, Base::G, Base::C]);
        let expected = Sequence::new(vec![0x18, 0x42], 4);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_try_from_byte_slice_for_sequence() {
        assert_eq!(
            Sequence::try_from(&b"ATCG"[..]),
            Ok(Sequence::from(vec![Base::A, Base::T, Base::C, Base::G]))
        );

        assert!(matches!(
            Sequence::try_from(&b"AT!G"[..]),
            Err(ParseError::InvalidBase(_))
        ));
    }

    #[test]
    fn test_from_sequence_for_sam_record_sequence() {
        use sam::record::{sequence::Base as SamBase, Sequence as SamSequence};

        let sequence = Sequence::new(vec![0x18, 0x42], 4); // ATGC

        let actual = SamSequence::from(&sequence);
        let expected = SamSequence::from(vec![SamBase::A, SamBase::T, SamBase::G, SamBase::C]);

        assert_eq!(actual, expected);
    }
}
