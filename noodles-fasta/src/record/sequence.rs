//! FASTA record sequence.

pub mod complement;

pub use self::complement::Complement;

use std::ops::Index;

use bytes::Bytes;
use noodles_core::{position::SequenceIndex, region::Interval};

/// A FASTA record sequence.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence(Bytes);

impl Sequence {
    /// Returns the length of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert_eq!(sequence.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether the sequence is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::default();
    /// assert!(sequence.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns a reference to a base at or slice of bases between the given index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_fasta::record::Sequence;
    ///
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    ///
    /// let start = Position::try_from(2)?;
    /// assert_eq!(sequence.get(start), Some(&b'C'));
    ///
    /// assert_eq!(sequence.get(start..), Some(&b"CGT"[..]));
    ///
    /// let end = Position::try_from(3)?;
    /// assert_eq!(sequence.get(start..=end), Some(&b"CG"[..]));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn get<I>(&self, index: I) -> Option<&I::Output>
    where
        I: SequenceIndex<u8>,
    {
        index.get(self.as_ref())
    }

    /// Returns a subset of the sequence within the given range.
    ///
    /// Unlike [`Self::get`], this returns the slice as a [`Sequence`].
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_fasta::record::Sequence;
    ///
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    ///
    /// let start = Position::try_from(2)?;
    /// let end = Position::try_from(3)?;
    /// let actual = sequence.slice(start..=end);
    ///
    /// let expected = Sequence::from(b"CG".to_vec());
    ///
    /// assert_eq!(actual, Some(expected));
    /// # Ok::<_, noodles_core::position::TryFromIntError>(())
    /// ```
    pub fn slice<I>(&self, interval: I) -> Option<Self>
    where
        I: Into<Interval>,
    {
        let interval = interval.into();

        let start = interval
            .start()
            .map(|position| usize::from(position) - 1)
            .unwrap_or(usize::MIN);

        let end = interval.end().map(usize::from).unwrap_or(self.len());

        if start <= end && end <= self.len() {
            let buf = self.0.slice(start..end);
            Some(Self::from(buf))
        } else {
            None
        }
    }

    /// Returns an iterator that complements the sequence.
    ///
    /// # Examples
    ///
    /// ## Complement a sequence
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let actual: Sequence = sequence.complement().collect::<Result<_, _>>()?;
    /// let expected = Sequence::from(b"TGCA".to_vec());
    /// assert_eq!(actual, expected);
    /// # Ok::<_, noodles_fasta::record::sequence::complement::ComplementError>(())
    /// ```
    ///
    /// ## Reverse complement a sequence
    ///
    /// ```
    /// use noodles_fasta::record::Sequence;
    /// let sequence = Sequence::from(b"ACGT".to_vec());
    /// let actual: Sequence = sequence.complement().rev().collect::<Result<_, _>>()?;
    /// let expected = sequence.clone();
    /// assert_eq!(actual, expected);
    /// # Ok::<_, noodles_fasta::record::sequence::complement::ComplementError>(())
    /// ```
    pub fn complement(&self) -> Complement<'_> {
        Complement::new(self.0.iter())
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl From<Vec<u8>> for Sequence {
    fn from(data: Vec<u8>) -> Self {
        Self(Bytes::from(data))
    }
}

impl From<Bytes> for Sequence {
    fn from(data: Bytes) -> Self {
        Self(data)
    }
}

impl FromIterator<u8> for Sequence {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = u8>,
    {
        let mut buf = Vec::new();
        buf.extend(iter);
        Self::from(buf)
    }
}

impl<I> Index<I> for Sequence
where
    I: SequenceIndex<u8>,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        index.index(self.as_ref())
    }
}
