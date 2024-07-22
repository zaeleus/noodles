/// An alignment record sequence.
pub trait Sequence {
    /// Returns whether there are any bases.
    fn is_empty(&self) -> bool;

    /// Returns the number of bases.
    fn len(&self) -> usize;

    /// Returns the base at the given index.
    fn get(&self, i: usize) -> Option<u8>;

    /// Splits the subsequence into two subsequences at the given index.
    fn split_at_checked(
        &self,
        mid: usize,
    ) -> Option<(Box<dyn Sequence + '_>, Box<dyn Sequence + '_>)>;

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

impl Sequence for Box<dyn Sequence + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }

    fn get(&self, i: usize) -> Option<u8> {
        (**self).get(i)
    }

    fn split_at_checked(
        &self,
        mid: usize,
    ) -> Option<(Box<dyn Sequence + '_>, Box<dyn Sequence + '_>)> {
        (**self).split_at_checked(mid)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        (**self).iter()
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

            fn get(&self, i: usize) -> Option<u8> {
                self.0.get(i).copied()
            }

            fn split_at_checked(
                &self,
                mid: usize,
            ) -> Option<(Box<dyn Sequence + '_>, Box<dyn Sequence + '_>)> {
                let (left, right) = self.0.split_at(mid);
                Some((Box::new(T(left.to_vec())), Box::new(T(right.to_vec()))))
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
