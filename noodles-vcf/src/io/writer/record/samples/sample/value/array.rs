use std::io::{self, Write};

use super::write_string;
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
                    write_string(writer, &s)?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::samples::sample::value::Array as ArrayBuf;

    #[test]
    fn test_write_array() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, array: &ArrayBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_array(buf, &array.into())?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let array = ArrayBuf::Integer(vec![Some(8)]);
        t(&mut buf, &array, b"8")?;

        let array = ArrayBuf::Integer(vec![Some(8), Some(13), None]);
        t(&mut buf, &array, b"8,13,.")?;

        let array = ArrayBuf::Float(vec![Some(0.333)]);
        t(&mut buf, &array, b"0.333")?;

        let array = ArrayBuf::Float(vec![Some(0.333), Some(0.667), None]);
        t(&mut buf, &array, b"0.333,0.667,.")?;

        let array = ArrayBuf::Character(vec![Some('n')]);
        t(&mut buf, &array, b"n")?;

        let array = ArrayBuf::Character(vec![Some('n'), Some('d'), None]);
        t(&mut buf, &array, b"n,d,.")?;

        let array = ArrayBuf::String(vec![Some(String::from("noodles"))]);
        t(&mut buf, &array, b"noodles")?;

        let array = ArrayBuf::String(vec![
            Some(String::from("noodles")),
            Some(String::from(";")),
            None,
        ]);
        t(&mut buf, &array, b"noodles,%3B,.")?;

        Ok(())
    }
}
