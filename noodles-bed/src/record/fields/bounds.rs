use std::ops::Range;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Bounds {
    pub(crate) reference_sequence_name_end: usize,
    pub(crate) feature_start_end: usize,
    pub(crate) feature_end_end: usize,
    pub(crate) other_fields_ends: Vec<usize>,
}

impl Bounds {
    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        0..self.reference_sequence_name_end
    }

    pub fn feature_start_range(&self) -> Range<usize> {
        self.reference_sequence_name_end..self.feature_start_end
    }

    pub fn feature_end_range(&self) -> Range<usize> {
        self.feature_start_end..self.feature_end_end
    }

    pub fn get(&self, i: usize) -> Option<Range<usize>> {
        let end = self.other_fields_ends.get(i).copied()?;

        let start = i
            .checked_sub(1)
            .and_then(|prev_i| self.other_fields_ends.get(prev_i).copied())
            .unwrap_or(self.feature_end_end);

        Some(start..end)
    }
}

impl Default for Bounds {
    fn default() -> Self {
        Self {
            reference_sequence_name_end: 3,
            feature_start_end: 4,
            feature_end_end: 5,
            other_fields_ends: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ranges() {
        let bounds = Bounds::default();
        assert_eq!(bounds.reference_sequence_name_range(), 0..3);
        assert_eq!(bounds.feature_start_range(), 3..4);
        assert_eq!(bounds.feature_end_range(), 4..5);
    }

    #[test]
    fn test_get() {
        let mut bounds = Bounds::default();

        assert!(bounds.get(0).is_none());

        bounds.other_fields_ends.push(6);
        bounds.other_fields_ends.push(7);

        assert_eq!(bounds.get(0), Some(5..6));
        assert_eq!(bounds.get(1), Some(6..7));
        assert!(bounds.get(2).is_none());
    }
}
