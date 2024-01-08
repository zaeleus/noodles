//! Alignment record CIGAR operation kind.

/// An alignment record CIGAR operation kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    /// An alignment match (`M`).
    Match,
    /// An insertion into the reference (`I`).
    Insertion,
    /// A deletion from the reference (`D`).
    Deletion,
    /// A skipped region from the reference (`N`).
    Skip,
    /// A soft clip (`S`).
    SoftClip,
    /// A hard clip (`H`).
    HardClip,
    /// Padding (`P`).
    Pad,
    /// A sequence match (`=`).
    SequenceMatch,
    /// A sequence mismatch (`X`).
    SequenceMismatch,
}

impl Kind {
    /// Returns whether the operation kind causes the alignment to consume the read.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::cigar::op::Kind;
    /// assert!(Kind::Match.consumes_read());
    /// assert!(Kind::Insertion.consumes_read());
    /// assert!(!Kind::Deletion.consumes_read());
    /// ```
    pub fn consumes_read(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Insertion
                | Self::SoftClip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }

    /// Returns whether the operation kind causes the alignment to consume the reference.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::alignment::record::cigar::op::Kind;
    /// assert!(Kind::Match.consumes_reference());
    /// assert!(!Kind::Insertion.consumes_reference());
    /// assert!(Kind::Deletion.consumes_reference());
    /// ```
    pub fn consumes_reference(&self) -> bool {
        matches!(
            self,
            Self::Match
                | Self::Deletion
                | Self::Skip
                | Self::SequenceMatch
                | Self::SequenceMismatch
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consumes_read() {
        assert!(Kind::Match.consumes_read());
        assert!(Kind::Insertion.consumes_read());
        assert!(!Kind::Deletion.consumes_read());
        assert!(!Kind::Skip.consumes_read());
        assert!(Kind::SoftClip.consumes_read());
        assert!(!Kind::HardClip.consumes_read());
        assert!(!Kind::Pad.consumes_read());
        assert!(Kind::SequenceMatch.consumes_read());
        assert!(Kind::SequenceMismatch.consumes_read());
    }

    #[test]
    fn test_consumes_reference() {
        assert!(Kind::Match.consumes_reference());
        assert!(!Kind::Insertion.consumes_reference());
        assert!(Kind::Deletion.consumes_reference());
        assert!(Kind::Skip.consumes_reference());
        assert!(!Kind::SoftClip.consumes_reference());
        assert!(!Kind::HardClip.consumes_reference());
        assert!(!Kind::Pad.consumes_reference());
        assert!(Kind::SequenceMatch.consumes_reference());
        assert!(Kind::SequenceMismatch.consumes_reference());
    }
}
