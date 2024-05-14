mod iter;

use std::fmt::{self, Write};

use noodles_sam as sam;

use self::iter::Iter;

/// A BAM record sequence.
#[derive(Eq, PartialEq)]
pub struct Sequence<'a> {
    src: &'a [u8],
    base_count: usize,
}

impl<'a> Sequence<'a> {
    pub(super) fn new(src: &'a [u8], base_count: usize) -> Self {
        Self { src, base_count }
    }

    /// Returns whether there are any bases.
    pub fn is_empty(&self) -> bool {
        self.src.is_empty()
    }

    /// Returns the number of bases in the sequence.
    ///
    /// This is _not_ the length of the buffer.
    pub fn len(&self) -> usize {
        self.base_count
    }

    /// Returns an iterator over the bases in the sequence.
    pub fn iter(&self) -> impl Iterator<Item = u8> + '_ {
        Iter::new(self.as_ref(), self.len())
    }
}

impl<'a> sam::alignment::record::Sequence for Sequence<'a> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.iter())
    }
}

impl<'a> AsRef<[u8]> for Sequence<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl<'a> fmt::Debug for Sequence<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        struct BasesFormat<'a>(&'a [u8], usize);

        impl<'a> fmt::Debug for BasesFormat<'a> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                use self::iter::decoded_bases;

                f.write_char('"')?;

                let BasesFormat(src, base_count) = *self;

                for b in src.iter().copied().flat_map(decoded_bases).take(base_count) {
                    f.write_char(char::from(b))?;
                }

                f.write_char('"')?;

                Ok(())
            }
        }

        f.debug_tuple("Sequence")
            .field(&BasesFormat(self.src, self.base_count))
            .finish()
    }
}

impl<'a> From<Sequence<'a>> for sam::alignment::record_buf::Sequence {
    fn from(sequence: Sequence<'a>) -> Self {
        Self::from(sequence.as_ref().to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_debug_fmt() {
        let sequence = Sequence::new(&[], 0);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("")"#);

        let sequence = Sequence::new(&[0x12, 0x40], 3);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("ACG")"#);

        let sequence = Sequence::new(&[0x12, 0x48], 4);
        assert_eq!(format!("{sequence:?}"), r#"Sequence("ACGT")"#);
    }
}
