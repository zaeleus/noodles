/// A SAM record sequence.
#[derive(Debug, Eq, PartialEq)]
pub struct Sequence<'a>(&'a [u8]);

impl<'a> Sequence<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Self(buf)
    }

    /// Returns whether the sequence is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the base at the given index.
    pub fn get(&self, i: usize) -> Option<u8> {
        self.0.get(i).copied()
    }

    /// Returns the number of bases in the sequence.
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl AsRef<[u8]> for Sequence<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl crate::alignment::record::Sequence for Sequence<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn len(&self) -> usize {
        self.len()
    }

    fn get(&self, i: usize) -> Option<u8> {
        self.get(i)
    }

    fn split_at_checked(
        &self,
        mid: usize,
    ) -> Option<(
        Box<dyn crate::alignment::record::Sequence + '_>,
        Box<dyn crate::alignment::record::Sequence + '_>,
    )> {
        if mid <= self.len() {
            let (left, right) = self.0.split_at(mid);
            Some((
                Box::new(Sequence::new(left)),
                Box::new(Sequence::new(right)),
            ))
        } else {
            None
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.as_ref().iter().copied())
    }
}

impl<'a> From<Sequence<'a>> for crate::alignment::record_buf::Sequence {
    fn from(sequence: Sequence<'a>) -> Self {
        Self::from(sequence.0.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crate_alignment_record_sequence_iter() {
        fn t(src: &[u8]) {
            let sequence: &dyn crate::alignment::record::Sequence = &Sequence::new(src);
            let actual: Vec<_> = sequence.iter().collect();
            assert_eq!(actual, src);
        }

        t(&[]);
        t(b"ACG");
        t(b"ACGT");
    }
}
