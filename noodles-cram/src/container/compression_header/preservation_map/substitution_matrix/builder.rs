use super::{Histogram, SubstitutionMatrix};
use crate::{
    record::{feature::substitution, Feature},
    Record,
};

#[derive(Debug, Default)]
pub struct Builder {
    histogram: Histogram,
}

impl Builder {
    pub fn update(&mut self, record: &Record) {
        for feature in record.features().iter() {
            if let Feature::Substitution { value, .. } = feature {
                match value {
                    substitution::Value::Bases(reference_base, read_base) => {
                        self.histogram.hit(*reference_base, *read_base);
                    }
                    substitution::Value::Code(_) => {
                        panic!("substitution matrix cannot be built from substitution codes");
                    }
                }
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
    use noodles_sam::alignment::record_buf::Sequence;

    use super::*;

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::{
            feature::substitution::{self, Base},
            Features,
        };

        // reference sequence = "ACAGGAATAANNNNNN"
        let bases: Sequence = Sequence::from(b"TCTGGCGTGT");

        let record = Record::builder()
            .set_alignment_start(Position::try_from(1)?)
            .set_read_length(bases.len())
            .set_bases(bases)
            .set_features(Features::from(vec![
                Feature::Substitution {
                    position: Position::try_from(1)?,
                    value: substitution::Value::Bases(Base::A, Base::T),
                },
                Feature::Substitution {
                    position: Position::try_from(3)?,
                    value: substitution::Value::Bases(Base::A, Base::T),
                },
                Feature::Substitution {
                    position: Position::try_from(6)?,
                    value: substitution::Value::Bases(Base::A, Base::C),
                },
                Feature::Substitution {
                    position: Position::try_from(7)?,
                    value: substitution::Value::Bases(Base::A, Base::G),
                },
                Feature::Substitution {
                    position: Position::try_from(9)?,
                    value: substitution::Value::Bases(Base::A, Base::G),
                },
                Feature::Substitution {
                    position: Position::try_from(10)?,
                    value: substitution::Value::Bases(Base::A, Base::T),
                },
            ]))
            .build();

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
