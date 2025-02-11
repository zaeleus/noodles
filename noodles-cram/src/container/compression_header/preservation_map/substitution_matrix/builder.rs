use super::{Histogram, SubstitutionMatrix};
use crate::{
    io::writer::{record::Feature, Record},
    record::feature::substitution::Base,
};

#[derive(Debug, Default)]
pub struct Builder {
    histogram: Histogram,
}

impl Builder {
    pub fn update(&mut self, record: &Record) {
        for feature in &record.features {
            if let Feature::Substitution {
                reference_base,
                read_base,
                ..
            } = feature
            {
                let reference_base = Base::try_from(*reference_base).unwrap();
                let read_base = Base::try_from(*read_base).unwrap();
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
    use noodles_core::Position;

    use super::*;

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        // reference sequence = "ACAGGAATAANNNNNN"
        let sequence = b"TCTGGCGTGT".to_vec();

        let record = Record {
            read_length: sequence.len(),
            alignment_start: Position::new(1),
            sequence,
            features: vec![
                Feature::Substitution {
                    position: Position::try_from(1)?,
                    reference_base: b'A',
                    read_base: b'T',
                },
                Feature::Substitution {
                    position: Position::try_from(3)?,
                    reference_base: b'A',
                    read_base: b'T',
                },
                Feature::Substitution {
                    position: Position::try_from(6)?,
                    reference_base: b'A',
                    read_base: b'C',
                },
                Feature::Substitution {
                    position: Position::try_from(7)?,
                    reference_base: b'A',
                    read_base: b'G',
                },
                Feature::Substitution {
                    position: Position::try_from(9)?,
                    reference_base: b'A',
                    read_base: b'G',
                },
                Feature::Substitution {
                    position: Position::try_from(10)?,
                    reference_base: b'A',
                    read_base: b'T',
                },
            ],
            ..Default::default()
        };

        let mut builder = Builder::default();
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

        Ok(())
    }
}
