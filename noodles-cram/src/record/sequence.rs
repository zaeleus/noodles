mod iter;

use noodles_core::Position;
use noodles_sam as sam;

use self::iter::Iter;
use super::Feature;
use crate::container::compression_header::preservation_map::SubstitutionMatrix;

pub(super) struct Sequence<'r, 'c: 'r> {
    reference_sequence: Option<&'c [u8]>,
    substitution_matrix: SubstitutionMatrix,
    features: &'r [Feature<'c>],
    alignment_start: Position,
    read_length: usize,
}

impl<'r, 'c: 'r> Sequence<'r, 'c> {
    pub(super) fn new(
        reference_sequence: Option<&'c [u8]>,
        substitution_matrix: SubstitutionMatrix,
        features: &'r [Feature<'c>],
        alignment_start: Position,
        read_length: usize,
    ) -> Self {
        Self {
            reference_sequence,
            substitution_matrix,
            features,
            alignment_start,
            read_length,
        }
    }
}

impl sam::alignment::record::Sequence for Sequence<'_, '_> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        self.read_length
    }

    fn get(&self, i: usize) -> Option<u8> {
        self.iter().nth(i)
    }

    fn split_at_checked(
        &self,
        mid: usize,
    ) -> Option<(
        Box<dyn sam::alignment::record::Sequence + '_>,
        Box<dyn sam::alignment::record::Sequence + '_>,
    )> {
        if mid > self.read_length {
            return None;
        }

        let mut left: Vec<u8> = self.iter().collect();
        let right = left.split_off(mid);
        Some((Box::new(OwnedBases(left)), Box::new(OwnedBases(right))))
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(Iter::new(
            self.reference_sequence,
            self.substitution_matrix.clone(),
            self.features,
            self.alignment_start,
            self.read_length,
        ))
    }
}

struct OwnedBases(Vec<u8>);

impl sam::alignment::record::Sequence for OwnedBases {
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
    ) -> Option<(
        Box<dyn sam::alignment::record::Sequence + '_>,
        Box<dyn sam::alignment::record::Sequence + '_>,
    )> {
        if mid <= self.0.len() {
            let (left, right) = self.0.split_at(mid);
            Some((
                Box::new(OwnedBases(left.to_vec())),
                Box::new(OwnedBases(right.to_vec())),
            ))
        } else {
            None
        }
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(self.0.iter().copied())
    }
}
