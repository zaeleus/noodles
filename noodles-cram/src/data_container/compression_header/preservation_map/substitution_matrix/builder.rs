use crate::{record::Feature, Record};

use super::{Base, Histogram, SubstitutionMatrix};

#[derive(Debug, Default)]
pub struct Builder {
    substitution_matrix: SubstitutionMatrix,
    histogram: Histogram,
}

impl Builder {
    pub fn update(&mut self, reference_sequence: &[u8], record: &Record) {
        for feature in record.features() {
            if let Feature::Substitution(pos, code) = feature {
                // FIXME: pos = 1-based, position = 0-based
                let reference_position = (pos - 1) as usize;
                let base = reference_sequence[reference_position] as char;
                let reference_base = Base::try_from(base).unwrap_or_default();
                let read_base = self.substitution_matrix.get(reference_base, *code);
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

        let mut builder = Builder::default();
        builder.update(reference_sequence, &record);
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
