/// Alignment record quality scores.
pub trait QualityScores {
    /// Returns whether there are any scores.
    fn is_empty(&self) -> bool;

    /// Returns the number of scores.
    fn len(&self) -> usize;

    /// Returns an iterator over scores.
    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_>;
}

impl<'a> IntoIterator for &'a dyn QualityScores {
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

        impl QualityScores for T {
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

        let quality_scores: &dyn QualityScores = &T(vec![45, 35, 43, 50]);

        assert_eq!(
            quality_scores.into_iter().collect::<Vec<_>>(),
            [45, 35, 43, 50]
        );
    }
}
