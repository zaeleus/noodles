use std::io;

/// Alignment record CIGAR operations.
pub trait Cigar {
    /// Returns whether there are any operations.
    fn is_empty(&self) -> bool;

    /// Returns the number of operations.
    fn len(&self) -> usize;

    /// Returns an iterator over operations.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(u8, usize)>> + '_>;
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

    #[test]
    fn test_into_iter() -> io::Result<()> {
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

        let sequence: &dyn Cigar = &T(vec![(b'M', 4)]);
        assert_eq!(
            sequence.into_iter().collect::<io::Result<Vec<_>>>()?,
            [(b'M', 4)]
        );

        Ok(())
    }
}
