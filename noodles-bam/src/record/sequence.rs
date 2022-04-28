//! BAM record sequence.

use std::str::FromStr;

use bytes::BytesMut;
use noodles_sam as sam;

/// A raw BAM record sequence buffer.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Sequence {
    pub(crate) buf: BytesMut,
    pub(crate) len: usize,
}

impl Sequence {
    /// Returns the number of bases in the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    ///
    /// let sequence = bam::record::Sequence::default();
    /// assert_eq!(sequence.len(), 0);
    ///
    /// let sequence: bam::record::Sequence = "ACGT".parse()?;
    /// assert_eq!(sequence.len(), 4);
    /// # Ok::<(), bam::record::sequence::ParseError>(())
    /// ```
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns whether the sequence is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    ///
    /// let sequence = bam::record::Sequence::default();
    /// assert!(sequence.is_empty());
    ///
    /// let sequence: bam::record::Sequence = "ACGT".parse()?;
    /// assert!(!sequence.is_empty());
    /// # Ok::<(), bam::record::sequence::ParseError>(())
    /// ```
    pub fn is_empty(&self) -> bool {
        self.buf.is_empty()
    }

    /// Removes all bases from the sequence.
    ///
    /// This does affect the capacity of the sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam as bam;
    ///
    /// let mut sequence: bam::record::Sequence = "ACGT".parse()?;
    /// assert!(!sequence.is_empty());
    ///
    /// sequence.clear();
    /// assert!(sequence.is_empty());
    /// # Ok::<(), bam::record::sequence::ParseError>(())
    /// ```
    pub fn clear(&mut self) {
        self.buf.clear();
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

/// An error returned when a raw alignment record sequence fails to parse.
pub type ParseError = sam::alignment::record::sequence::ParseError;

impl FromStr for Sequence {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sam_sequence = s.parse()?;
        Ok(Self::from(&sam_sequence))
    }
}

impl From<&sam::alignment::record::Sequence> for Sequence {
    fn from(sequence: &sam::alignment::record::Sequence) -> Self {
        use crate::writer::record::put_sequence;

        let len = sequence.len();
        let mut buf = BytesMut::with_capacity((len + 1) / 2);
        put_sequence(&mut buf, sequence);

        Sequence { buf, len }
    }
}

impl From<Sequence> for sam::alignment::record::Sequence {
    fn from(bam_sequence: Sequence) -> Self {
        use crate::reader::alignment_record::get_sequence;

        let mut data = &bam_sequence.buf[..];
        let mut sam_sequence = sam::alignment::record::Sequence::default();
        // FIXME
        get_sequence(&mut data, &mut sam_sequence, bam_sequence.len()).unwrap();

        sam_sequence
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sam_alignment_record_sequence_for_sequence(
    ) -> Result<(), sam::alignment::record::sequence::ParseError> {
        let sam_sequence = "ACGT".parse()?;
        let actual = Sequence::from(&sam_sequence);

        let expected = Sequence {
            buf: BytesMut::from(&[0x12, 0x48][..]),
            len: 4,
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_from_sequence_for_sam_alignment_record_sequence(
    ) -> Result<(), sam::alignment::record::sequence::ParseError> {
        let bam_sequence = Sequence {
            buf: BytesMut::from(&[0x12, 0x48][..]),
            len: 4,
        };
        let actual = sam::alignment::record::Sequence::from(bam_sequence);

        let expected = "ACGT".parse()?;

        assert_eq!(actual, expected);

        Ok(())
    }
}
