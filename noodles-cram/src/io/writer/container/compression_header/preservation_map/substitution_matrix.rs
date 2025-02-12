mod frequencies;

use self::frequencies::Frequencies;
use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    io::writer::{record::Feature, Record},
};

pub(super) fn build_substitution_matrix(records: &[Record]) -> SubstitutionMatrix {
    let mut frequencies = Frequencies::default();

    for record in records {
        for feature in &record.features {
            if let Feature::Substitution {
                reference_base,
                read_base,
                ..
            } = feature
            {
                frequencies.hit(*reference_base, *read_base);
            }
        }
    }

    frequencies.into()
}
