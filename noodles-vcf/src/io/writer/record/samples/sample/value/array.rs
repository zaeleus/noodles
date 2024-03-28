use std::io::{self, Write};

use crate::{io::writer::record::MISSING, variant::record::samples::series::value::Array};

pub(super) fn write_array<W>(writer: &mut W, array: &Array) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b",";

    match array {
        Array::Integer(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = result? {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::Float(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = result? {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::Character(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(c) = result? {
                    write!(writer, "{c}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
        Array::String(values) => {
            for (i, result) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(s) = result? {
                    writer.write_all(s.as_bytes())?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
    }

    Ok(())
}
