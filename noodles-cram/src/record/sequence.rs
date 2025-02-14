mod iter;

use noodles_core::Position;
use noodles_fasta as fasta;
use noodles_sam as sam;

use self::iter::Iter;
use super::Feature;
use crate::container::compression_header::preservation_map::SubstitutionMatrix;

pub(super) struct Sequence<'r, 'c: 'r> {
    reference_sequence: Option<fasta::record::Sequence>,
    substitution_matrix: SubstitutionMatrix,
    features: &'r [Feature<'c>],
    alignment_start: Position,
    read_length: usize,
}

impl<'r, 'c: 'r> Sequence<'r, 'c> {
    pub(super) fn new(
        reference_sequence: Option<fasta::record::Sequence>,
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
        _mid: usize,
    ) -> Option<(
        Box<dyn sam::alignment::record::Sequence + '_>,
        Box<dyn sam::alignment::record::Sequence + '_>,
    )> {
        todo!()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = u8> + '_> {
        Box::new(Iter::new(
            self.reference_sequence.as_ref(),
            self.substitution_matrix.clone(),
            self.features,
            self.alignment_start,
            self.read_length,
        ))
    }
}
