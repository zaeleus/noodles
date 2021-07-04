//! BAM CIGAR and operations.

pub mod op;
mod ops;

pub use self::{op::Op, ops::Ops};

use std::{convert::TryFrom, fmt, io, ops::Deref};

use noodles_sam::{self as sam, record::cigar::op::Kind};

/// BAM record CIGAR.
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    /// Creates a CIGAR by wrapping raw CIGAR data.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::Cigar;
    ///
    /// let data = [
    ///     0x40, 0x02, 0x00, 0x00, // 36M
    ///     0x84, 0x00, 0x00, 0x00, // 8S
    /// ];
    ///
    /// let cigar = Cigar::new(&data);
    ///
    /// assert_eq!(*cigar, data);
    /// ```
    pub fn new(bytes: &[u8]) -> Cigar<'_> {
        Cigar(bytes)
    }

    /// Returns a iterator over the operations in the CIGAR.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{cigar::Op, Cigar};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// let data = [
    ///     0x40, 0x02, 0x00, 0x00, // 36M
    ///     0x84, 0x00, 0x00, 0x00, // 8S
    /// ];
    ///
    /// let cigar = Cigar::new(&data);
    ///
    /// let mut ops = cigar.ops();
    ///
    /// assert_eq!(ops.next().transpose()?, Some(Op::new(Kind::Match, 36)?));
    /// assert_eq!(ops.next().transpose()?, Some(Op::new(Kind::SoftClip, 8)?));
    /// assert_eq!(ops.next().transpose()?, None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn ops(&self) -> Ops<'_> {
        Ops::new(self.0)
    }

    /// Calculates the alignment span over the reference sequence.
    ///
    /// This sums the lengths of the CIGAR operations that consume the reference sequence, i.e.,
    /// alignment matches (`M`), deletions from the reference (`D`), skipped reference regions
    /// (`S`), sequence matches (`=`), and sequence mismatches (`X`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bam::record::{cigar::Op, Cigar};
    /// use noodles_sam::record::cigar::op::Kind;
    ///
    /// let data = [
    ///     0x40, 0x02, 0x00, 0x00, // 36M
    ///     0x43, 0x00, 0x00, 0x00, // 4D
    ///     0x84, 0x00, 0x00, 0x00, // 8S
    /// ];
    ///
    /// let cigar = Cigar::new(&data);
    ///
    /// assert_eq!(cigar.reference_len()?, 40);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn reference_len(&self) -> io::Result<u32> {
        let mut len = 0;

        for result in self.ops() {
            let op = result?;

            match op.kind() {
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SeqMatch | Kind::SeqMismatch => {
                    len += op.len();
                }
                _ => {}
            }
        }

        Ok(len)
    }
}

impl<'a> fmt::Debug for Cigar<'a> {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt.debug_list().entries(self.ops()).finish()
    }
}

impl<'a> fmt::Display for Cigar<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for result in self.ops() {
            let op = result.map_err(|_| fmt::Error)?;
            write!(f, "{}", op)?;
        }

        Ok(())
    }
}

impl<'a> Deref for Cigar<'a> {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Cigar<'a>> for sam::record::Cigar {
    type Error = io::Error;

    fn try_from(cigar: Cigar<'_>) -> Result<Self, Self::Error> {
        let mut ops = Vec::new();

        for result in cigar.ops() {
            let op = result?;
            ops.push(sam::record::cigar::Op::new(op.kind(), op.len()));
        }

        Ok(Self::from(ops))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bytes() -> Result<(), Box<dyn std::error::Error>> {
        let bytes = [0x40, 0x02, 0x00, 0x00, 0x62, 0x03, 0x00, 0x00];
        let cigar = Cigar::new(&bytes);

        let mut ops = cigar.ops();

        assert_eq!(ops.next().transpose()?, Some(Op::try_from(0x240)?));
        assert_eq!(ops.next().transpose()?, Some(Op::try_from(0x362)?));
        assert_eq!(ops.next().transpose()?, None);

        Ok(())
    }

    #[test]
    fn test_try_from_cigar_for_sam_record_cigar() -> io::Result<()> {
        use sam::record::cigar::{op, Op};

        let bytes = [0x40, 0x02, 0x00, 0x00, 0x62, 0x03, 0x00, 0x00];
        let cigar = Cigar::new(&bytes);

        let actual = sam::record::Cigar::try_from(cigar)?;
        let expected = sam::record::Cigar::from(vec![
            Op::new(op::Kind::Match, 36),
            Op::new(op::Kind::Deletion, 54),
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }
}
