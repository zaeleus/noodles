use std::io::{self, Write};

use crate::variant::record::samples::series::value::{genotype::Phasing, Genotype};

pub(super) fn write_genotype<W>(writer: &mut W, genotype: &dyn Genotype) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = b'.';

    for (i, result) in genotype.iter().enumerate() {
        let (position, phasing) = result?;

        if i > 0 {
            write_phasing(writer, phasing)?;
        }

        if let Some(n) = position {
            write!(writer, "{n}")?;
        } else {
            writer.write_all(&[MISSING])?;
        }
    }

    Ok(())
}

fn write_phasing<W>(writer: &mut W, phasing: Phasing) -> io::Result<()>
where
    W: Write,
{
    const PHASED: u8 = b'/';
    const UNPHASED: u8 = b'|';

    match phasing {
        Phasing::Phased => writer.write_all(&[PHASED]),
        Phasing::Unphased => writer.write_all(&[UNPHASED]),
    }
}
