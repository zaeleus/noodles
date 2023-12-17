use std::io;

/// Alignment record CIGAR operations.
pub trait Cigar {
    /// Returns whether there are any operations.
    fn is_empty(&self) -> bool;

    /// Returns the number of operations.
    fn len(&self) -> usize;

    /// Returns an iterator over operations.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(u8, usize)>> + '_>;

    /// Calculates the alignment span over the reference sequence.
    fn alignment_span(&self) -> io::Result<usize> {
        fn consumes_reference(op: u8) -> bool {
            matches!(op, b'M' | b'D' | b'N' | b'=' | b'X')
        }

        let mut span = 0;

        for result in self.iter() {
            let (op, len) = result?;

            if consumes_reference(op) {
                span += len;
            }
        }

        Ok(span)
    }
}

impl<'a> IntoIterator for &'a dyn Cigar {
    type Item = io::Result<(u8, usize)>;
    type IntoIter = Box<dyn Iterator<Item = Self::Item> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct T(Vec<(u8, usize)>);

    impl Cigar for T {
        fn is_empty(&self) -> bool {
            self.0.is_empty()
        }

        fn len(&self) -> usize {
            self.0.len()
        }

        fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(u8, usize)>> + '_> {
            Box::new(self.0.iter().copied().map(Ok))
        }
    }

    #[test]
    fn test_into_iter() -> io::Result<()> {
        let cigar: &dyn Cigar = &T(vec![(b'M', 4)]);

        assert_eq!(
            cigar.into_iter().collect::<io::Result<Vec<_>>>()?,
            [(b'M', 4)]
        );

        Ok(())
    }

    #[test]
    fn test_alignment_span() -> io::Result<()> {
        let cigar: &dyn Cigar = &T(vec![(b'M', 36), (b'D', 4), (b'S', 8)]);
        assert_eq!(cigar.alignment_span()?, 40);
        Ok(())
    }
}
