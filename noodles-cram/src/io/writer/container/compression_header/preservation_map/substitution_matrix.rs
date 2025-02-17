mod frequencies;

use std::io::{self, Write};

use self::frequencies::Frequencies;
use crate::{
    container::compression_header::preservation_map::SubstitutionMatrix,
    io::writer::{record::Feature, Record},
};

pub(super) fn write_substitution_matrix<W>(
    writer: &mut W,
    substitution_matrix: &SubstitutionMatrix,
) -> io::Result<()>
where
    W: Write,
{
    let buf = <[u8; 5]>::from(substitution_matrix);
    writer.write_all(&buf)
}

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
