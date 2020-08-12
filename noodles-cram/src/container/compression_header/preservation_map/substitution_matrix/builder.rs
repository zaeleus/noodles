use std::convert::TryFrom;

use crate::{record::Feature, Record};

use super::{Base, Histogram, SubstitutionMatrix};

#[derive(Debug)]
pub struct Builder<'a> {
    reference_sequence: &'a [u8],
    substitution_matrix: SubstitutionMatrix,
    histogram: Histogram,
}

impl<'a> Builder<'a> {
    pub fn new(reference_sequence: &'a [u8]) -> Self {
        Self {
            reference_sequence,
            substitution_matrix: SubstitutionMatrix::default(),
            histogram: Histogram::default(),
        }
    }

    pub fn update(&mut self, record: &Record) {
        let substitution_matrix = SubstitutionMatrix::default();

        for feature in &record.features {
            if let Feature::Substitution(pos, code) = feature {
                // FIXME: pos = 1-based, position = 0-based
                let reference_position = (pos - 1) as usize;
                let base = self.reference_sequence[reference_position] as char;
                let reference_base = Base::try_from(base).unwrap_or_default();
                let read_base = substitution_matrix.get(reference_base, *code);
                self.histogram.hit(reference_base, read_base);
            }
        }
    }

    pub fn build(self) -> SubstitutionMatrix {
        SubstitutionMatrix::from(self.histogram)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build() {
        let reference_sequence = b"ACAGGAATAANNNNNN";

        let mut record = Record::default();
        record.add_feature(Feature::Substitution(1, 2)); // A => T
        record.add_feature(Feature::Substitution(3, 2)); // A => T
        record.add_feature(Feature::Substitution(6, 0)); // A => C
        record.add_feature(Feature::Substitution(7, 1)); // A => G
        record.add_feature(Feature::Substitution(9, 1)); // A => G
        record.add_feature(Feature::Substitution(10, 2)); // A => T

        let mut builder = Builder::new(reference_sequence);
        builder.update(&record);
        let matrix = builder.build();

        assert_eq!(
            matrix.substitutions,
            [
                [Base::T, Base::G, Base::C, Base::N],
                [Base::A, Base::G, Base::T, Base::N],
                [Base::A, Base::C, Base::T, Base::N],
                [Base::A, Base::C, Base::G, Base::N],
                [Base::A, Base::C, Base::G, Base::T],
            ]
        );
    }
}
