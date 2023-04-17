use std::io::{self, Write};

use super::MISSING;
use crate::record::{info::field::Value, Info};

pub(super) fn write_info<W>(writer: &mut W, info: &Info) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b";";
    const SEPARATOR: &[u8] = b"=";

    if info.is_empty() {
        writer.write_all(MISSING)?;
    } else {
        for (i, (key, value)) in info.as_ref().iter().enumerate() {
            if i > 0 {
                writer.write_all(DELIMITER)?;
            }

            writer.write_all(key.as_ref().as_bytes())?;

            match value {
                Some(Value::Flag) => {}
                Some(v) => {
                    writer.write_all(SEPARATOR)?;
                    write_value(writer, v)?;
                }
                None => {
                    writer.write_all(SEPARATOR)?;
                    writer.write_all(MISSING)?;
                }
            }
        }
    }

    Ok(())
}

fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    const DELIMITER: &[u8] = b",";

    match value {
        Value::Integer(n) => write!(writer, "{n}"),
        Value::Float(n) => write!(writer, "{n}"),
        Value::Flag => Ok(()),
        Value::Character(c) => write!(writer, "{c}"),
        Value::String(s) => writer.write_all(s.as_bytes()),
        Value::IntegerArray(values) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = v {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::FloatArray(values) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(n) = v {
                    write!(writer, "{n}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::CharacterArray(values) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(c) = v {
                    write!(writer, "{c}")?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
        Value::StringArray(values) => {
            for (i, v) in values.iter().enumerate() {
                if i > 0 {
                    writer.write_all(DELIMITER)?;
                }

                if let Some(s) = v {
                    writer.write_all(s.as_bytes())?;
                } else {
                    writer.write_all(MISSING)?;
                }
            }

            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_info() -> io::Result<()> {
        use crate::record::info::field::key;

        fn t(buf: &mut Vec<u8>, info: &Info, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info(buf, info)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let info = Info::default();
        t(&mut buf, &info, b".")?;

        let info = [(key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2)))]
            .into_iter()
            .collect();
        t(&mut buf, &info, b"NS=2")?;

        let info = [
            (key::SAMPLES_WITH_DATA_COUNT, Some(Value::Integer(2))),
            (key::IS_IN_DB_SNP, Some(Value::Flag)),
        ]
        .into_iter()
        .collect();

        t(&mut buf, &info, b"NS=2;DB")?;

        Ok(())
    }
}
