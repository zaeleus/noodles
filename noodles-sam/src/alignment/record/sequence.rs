/// An alignment record sequence.
pub trait Sequence {
    /// Returns whether there are any bases.
    fn is_empty(&self) -> bool;

    /// Returns the number of bases.
    fn len(&self) -> usize;

    /// Returns an iterator over bases.
    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_>;
}

impl<'a> IntoIterator for &'a dyn Sequence {
    type Item = u8;
    type IntoIter = Box<dyn Iterator<Item = Self::Item> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_into_iter() {
        struct T(Vec<u8>);

        impl Sequence for T {
            fn is_empty(&self) -> bool {
                self.0.is_empty()
            }

            fn len(&self) -> usize {
                self.0.len()
            }

            fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
                Box::new(self.0.iter().copied())
            }
        }

        let sequence: &dyn Sequence = &T(vec![b'N', b'D', b'L', b'S']);

        assert_eq!(
            sequence.into_iter().collect::<Vec<_>>(),
            [b'N', b'D', b'L', b'S']
        );
    }
}
