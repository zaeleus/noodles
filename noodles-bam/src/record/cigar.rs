use std::{fmt, io, mem};

use noodles_sam::{self as sam, alignment::record::cigar::Op};

const CHUNK_SIZE: usize = mem::size_of::<u32>();

/// BAM record CIGAR operations.
#[derive(Eq, PartialEq)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any CIGAR operations.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the number of CIGAR operations.
    ///
    /// This is _not_ the length of the buffer.
    pub fn len(&self) -> usize {
        self.0.len() / CHUNK_SIZE
    }

    /// Returns an iterator over CIGAR operations.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<Op>> + '_ {
        use crate::record::codec::decoder::cigar::decode_op;

        self.0.chunks(CHUNK_SIZE).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let n = u32::from_le_bytes(buf);
            decode_op(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}

impl sam::alignment::record::Cigar for Cigar<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Op>> + '_> {
        Box::new(self.iter())
    }
}

impl AsRef<[u8]> for Cigar<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl fmt::Debug for Cigar<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<'a> TryFrom<Cigar<'a>> for sam::alignment::record_buf::Cigar {
    type Error = io::Error;

    fn try_from(bam_cigar: Cigar<'a>) -> Result<Self, Self::Error> {
        bam_cigar.iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let src = &[][..];
        let cigar = Cigar::new(src);
        assert!(cigar.is_empty());

        let src = &[0x40, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        assert!(!cigar.is_empty());

        let src = &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        assert!(!cigar.is_empty());
    }

    #[test]
    fn test_len() {
        let src = &[][..];
        let cigar = Cigar::new(src);
        assert_eq!(cigar.len(), 0);

        let src = &[0x40, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        assert_eq!(cigar.len(), 1);

        let src = &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        assert_eq!(cigar.len(), 2);
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        use sam::alignment::record::cigar::op::Kind;

        let src = &[][..];
        let cigar = Cigar::new(src);
        assert!(cigar.iter().next().is_none());

        let src = &[0x40, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        let actual: Vec<_> = cigar.iter().collect::<io::Result<_>>()?;
        let expected = [Op::new(Kind::Match, 4)];
        assert_eq!(actual, expected);

        let src = &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        let actual: Vec<_> = cigar.iter().collect::<io::Result<_>>()?;
        let expected = [Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)];
        assert_eq!(actual, expected);

        let src = &[0x40, 0x00, 0x00][..];
        let cigar = Cigar::new(src);
        let mut iter = cigar.iter();
        assert!(matches!(
            iter.next(),
            Some(Err(e)) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
