use std::{
    convert::TryFrom,
    io::{self, Read},
};

use noodles_vcf as vcf;

use crate::{header::StringMap, reader::value::read_value, record::Value};

pub fn read_info<R>(
    reader: &mut R,
    infos: &vcf::header::Infos,
    string_map: &StringMap,
    len: usize,
) -> io::Result<vcf::record::Info>
where
    R: Read,
{
    use vcf::record::info::Field;

    let mut fields = Vec::with_capacity(len);

    for _ in 0..len {
        let key = read_info_key(reader, string_map)?;

        let info = infos.get(&key).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing header INFO record for {}", key),
            )
        })?;

        let value = read_info_value(reader, &info)?;

        let field = Field::new(key, value);
        fields.push(field);
    }

    vcf::record::Info::try_from(fields).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn read_info_key<R>(
    reader: &mut R,
    string_map: &StringMap,
) -> io::Result<vcf::record::info::field::Key>
where
    R: Read,
{
    use crate::record::value::Int8;

    match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => usize::try_from(i)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|j| {
                string_map.get_index(j).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid string map index: {}", j),
                    )
                })
            })
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected i8, got {:?}", v),
            ));
        }
    }
}

fn read_info_value<R>(
    reader: &mut R,
    info: &vcf::header::Info,
) -> io::Result<vcf::record::info::field::Value>
where
    R: Read,
{
    use vcf::{header::info::Type, record::info};

    use crate::record::value::{Float, Int16, Int32, Int8};

    let value = match info.ty() {
        Type::Integer => match read_value(reader)? {
            Some(Value::Int8(Some(Int8::Value(n)))) => info::field::Value::Integer(i32::from(n)),
            Some(Value::Int8Array(values)) => {
                info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
            }
            Some(Value::Int16(Some(Int16::Value(n)))) => info::field::Value::Integer(i32::from(n)),
            Some(Value::Int16Array(values)) => {
                info::field::Value::IntegerArray(values.into_iter().map(i32::from).collect())
            }
            Some(Value::Int32(Some(Int32::Value(n)))) => info::field::Value::Integer(n),
            Some(Value::Int32Array(values)) => info::field::Value::IntegerArray(values),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Integer, v),
                ));
            }
        },
        Type::Flag => match read_value(reader)? {
            Some(Value::Int8(Some(Int8::Value(1)))) | None => info::field::Value::Flag,
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Flag, v),
                ));
            }
        },
        Type::Float => match read_value(reader)? {
            Some(Value::Float(Some(Float::Value(n)))) => info::field::Value::Float(n),
            Some(Value::FloatArray(values)) => info::field::Value::FloatArray(values),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Float, v),
                ))
            }
        },
        Type::Character => match read_value(reader)? {
            Some(Value::String(Some(s))) => s
                .chars()
                .next()
                .map(info::field::Value::Character)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "INFO character value missing")
                })?,
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::Character, v),
                ))
            }
        },
        Type::String => match read_value(reader)? {
            Some(Value::String(Some(s))) => info::field::Value::String(s),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("type mismatch: expected {}, got {:?}", Type::String, v),
                ))
            }
        },
    };

    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_info_value() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::{
            header::{info::Type, Info, Number},
            record,
        };

        let data = [
            0x17, // Type::String(1)
            0x6e, // "n"
        ];
        let mut reader = &data[..];

        let key = record::info::field::Key::Other(
            String::from("CHAR"),
            Number::Count(1),
            Type::Character,
            String::default(),
        );
        let info = Info::from(key);

        let actual = read_info_value(&mut reader, &info)?;
        let expected = record::info::field::Value::Character('n');
        assert_eq!(actual, expected);

        Ok(())
    }
}
