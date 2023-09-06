mod key;
mod values;

use std::io;

use noodles_vcf::{
    self as vcf,
    record::{genotypes::Keys, Genotypes},
};

use self::{
    key::read_key,
    values::{read_genotype_values, read_values},
};
use crate::header::string_maps::StringStringMap;

pub fn read_genotypes(
    src: &mut &[u8],
    formats: &vcf::header::Formats,
    string_map: &StringStringMap,
    sample_count: usize,
    format_count: usize,
) -> io::Result<Genotypes> {
    use vcf::record::genotypes::keys::key;

    let mut keys = Vec::with_capacity(format_count);
    let mut values = vec![Vec::new(); sample_count];

    for _ in 0..format_count {
        let key = read_key(src, formats, string_map)?;

        let vs = if key == &key::GENOTYPE {
            read_genotype_values(src, sample_count)?
        } else {
            read_values(src, sample_count)?
        };

        keys.push(key.clone());

        for (sample, value) in values.iter_mut().zip(vs) {
            sample.push(value);
        }
    }

    let keys = Keys::try_from(keys).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Genotypes::new(keys, values))
}
