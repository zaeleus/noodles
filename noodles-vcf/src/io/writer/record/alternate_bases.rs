use std::io::{self, Write};

use super::MISSING;
use crate::variant::record_buf::AlternateBases;

pub(super) fn write_alternate_bases<W>(
    writer: &mut W,
    alternate_bases: &AlternateBases,
) -> io::Result<()>
where
    W: Write,
{
    if alternate_bases.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, allele) in alternate_bases.as_ref().iter().enumerate() {
            if i > 0 {
                write!(writer, ",")?;
            }

            write!(writer, "{allele}")?;
        }
    }

    Ok(())
}
