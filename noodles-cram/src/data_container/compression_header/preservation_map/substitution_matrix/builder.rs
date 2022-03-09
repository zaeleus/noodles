use noodles_fasta as fasta;
use noodles_sam::AlignmentRecord;

use super::{Base, Histogram, SubstitutionMatrix};
use crate::{record::Feature, Record};

#[derive(Debug, Default)]
pub struct Builder {
    histogram: Histogram,
}

impl Builder {
    pub fn update(&mut self, reference_sequence: &fasta::record::Sequence, record: &Record) {
        let alignment_start = match record.alignment_start() {
            Some(position) => position,
            None => return,
        };

        let read_bases = record.bases();

        for ((reference_position, mut read_position), feature) in
            record.features().with_positions(alignment_start)
        {
            read_position -= 1;

            if let Feature::Substitution(..) = feature {
                let base = reference_sequence.get(reference_position).copied().unwrap();
                let reference_base = Base::try_from(base).unwrap_or_default();

                let base = read_bases[read_position];
                let read_base = Base::try_from(base).unwrap_or_default();

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
    use noodles_sam as sam;

    use super::*;

    #[test]
    fn test_build() -> Result<(), sam::record::position::TryFromIntError> {
        let reference_sequence = fasta::record::Sequence::from(b"ACAGGAATAANNNNNN".to_vec());
        let bases = b"TCTGGCGTGT";

        let record = Record::builder()
            .set_alignment_start(sam::record::Position::try_from(1)?)
            .set_read_length(bases.len())
            .set_bases(bases.to_vec())
            .add_feature(Feature::Substitution(1, 2)) // A => T
            .add_feature(Feature::Substitution(3, 2)) // A => T
            .add_feature(Feature::Substitution(6, 0)) // A => C
            .add_feature(Feature::Substitution(7, 1)) // A => G
            .add_feature(Feature::Substitution(9, 1)) // A => G
            .add_feature(Feature::Substitution(10, 2)) // A => T
            .build();

        let mut builder = Builder::default();
        builder.update(&reference_sequence, &record);
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

        Ok(())
    }
}
