mod value;

use std::io::{self, Write};

use self::value::write_value;
use crate::{io::writer::record::MISSING, variant::record::samples::Sample, Header};

pub(super) fn write_sample<W, S>(writer: &mut W, header: &Header, sample: S) -> io::Result<()>
where
    W: Write,
    S: Sample,
{
    const DELIMITER: &[u8] = b":";

    for (i, result) in sample.iter(header).enumerate() {
        let (_, value) = result?;

        if i > 0 {
            writer.write_all(DELIMITER)?;
        }

        match value {
            Some(v) => write_value(writer, header, &v)?,
            None => writer.write_all(MISSING)?,
        }
    }

    Ok(())
}
