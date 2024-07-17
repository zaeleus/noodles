use std::ops::Range;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Bounds<const N: usize> {
    pub(crate) standard_fields_ends: [usize; N],
    pub(crate) other_fields_ends: Vec<usize>,
}

impl<const N: usize> Bounds<N> {
    pub fn reference_sequence_name_range(&self) -> Range<usize> {
        0..self.standard_fields_ends[0]
    }

    pub fn feature_start_range(&self) -> Range<usize> {
        self.standard_fields_ends[0]..self.standard_fields_ends[1]
    }

    pub fn feature_end_range(&self) -> Range<usize> {
        self.standard_fields_ends[1]..self.standard_fields_ends[2]
    }

    pub fn get(&self, i: usize) -> Option<Range<usize>> {
        let end = self.other_fields_ends.get(i).copied()?;

        let start = i
            .checked_sub(1)
            .and_then(|prev_i| self.other_fields_ends.get(prev_i).copied())
            .unwrap_or(self.standard_fields_ends[N - 1]);

        Some(start..end)
    }
}

impl Bounds<4> {
    pub fn name_range(&self) -> Range<usize> {
        self.standard_fields_ends[2]..self.standard_fields_ends[3]
    }
}

impl Bounds<5> {
    pub fn name_range(&self) -> Range<usize> {
        self.standard_fields_ends[2]..self.standard_fields_ends[3]
    }

    pub fn score_range(&self) -> Range<usize> {
        self.standard_fields_ends[3]..self.standard_fields_ends[4]
    }
}

impl Default for Bounds<3> {
    fn default() -> Self {
        Self {
            standard_fields_ends: [3, 4, 5],
            other_fields_ends: Vec::new(),
        }
    }
}

impl Default for Bounds<4> {
    fn default() -> Self {
        Self {
            standard_fields_ends: [3, 4, 5, 6],
            other_fields_ends: Vec::new(),
        }
    }
}

impl Default for Bounds<5> {
    fn default() -> Self {
        Self {
            standard_fields_ends: [3, 4, 5, 6, 7],
            other_fields_ends: Vec::new(),
        }
    }
}

impl Default for Bounds<6> {
    fn default() -> Self {
        Self {
            standard_fields_ends: [3, 4, 5, 6, 7, 8],
            other_fields_ends: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ranges() {
        let bounds = Bounds::<3>::default();
        assert_eq!(bounds.reference_sequence_name_range(), 0..3);
        assert_eq!(bounds.feature_start_range(), 3..4);
        assert_eq!(bounds.feature_end_range(), 4..5);
    }

    #[test]
    fn test_get() {
        let mut bounds = Bounds::<3>::default();

        assert!(bounds.get(0).is_none());

        bounds.other_fields_ends.push(6);
        bounds.other_fields_ends.push(7);

        assert_eq!(bounds.get(0), Some(5..6));
        assert_eq!(bounds.get(1), Some(6..7));
        assert!(bounds.get(2).is_none());
    }
}
