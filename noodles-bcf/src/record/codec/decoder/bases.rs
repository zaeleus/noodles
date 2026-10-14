use std::{io, num::NonZero};

use noodles_vcf::variant::record_buf::AlternateBases;

use super::read_value;
use crate::record::codec::Value;

pub(crate) fn read_ref_alt(
    src: &mut &[u8],
    len: NonZero<usize>,
) -> io::Result<(String, AlternateBases)> {
    let mut alleles = Vec::with_capacity(len.get());

    for _ in 0..len.get() {
        match read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
            Some(Value::String(Some(s))) => alleles.push(s),
            Some(Value::String(None)) => alleles.push("."),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid ref_alt: expected string, got {v:?}"),
                ));
            }
        }
    }

    // SAFETY: `len > 0`.
    let (raw_reference_bases, raw_alternate_bases) = alleles.split_at(1);

    let reference_bases = raw_reference_bases
        .first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing reference bases"))
        .map(|&s| s.into())?;

    let alternate_bases = raw_alternate_bases
        .iter()
        .map(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        })
        .collect::<Result<_, _>>()?;

    Ok((reference_bases, alternate_bases))
}
