use std::{
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_vcf::record::genotypes::sample::Value;

use crate::{
    lazy::record::value::{Float, Int16, Int32, Int8, Type},
    record::codec::decoder::value::read_type,
};

pub(super) fn read_values(src: &mut &[u8], sample_count: usize) -> io::Result<Vec<Option<Value>>> {
    match read_type(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
        Some(Type::Int8(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_int8_values(src, sample_count),
            _ => read_int8_array_values(src, sample_count, len),
        },
        Some(Type::Int16(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_int16_values(src, sample_count),
            _ => read_int16_array_values(src, sample_count, len),
        },
        Some(Type::Int32(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_int32_values(src, sample_count),
            _ => read_int32_array_values(src, sample_count, len),
        },
        Some(Type::Float(len)) => match len {
            0 => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid number: {len}"),
            )),
            1 => read_float_values(src, sample_count),
            _ => read_float_array_values(src, sample_count, len),
        },
        Some(Type::String(len)) => read_string_values(src, sample_count, len),
        ty => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unhandled type: {ty:?}"),
        )),
    }
}

fn read_int8_values(src: &mut &[u8], sample_count: usize) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = src.read_i8().map(Int8::from)?;

        match value {
            Int8::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int8::Missing => values.push(None),
            _ => todo!("unhandled i8 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int8_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        src.read_i8_into(&mut buf)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int8::from)
            .filter_map(|value| match value {
                Int8::Value(n) => Some(Some(i32::from(n))),
                Int8::Missing => Some(None),
                Int8::EndOfVector => None,
                _ => todo!("unhandled i8 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_int16_values(src: &mut &[u8], sample_count: usize) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = src.read_i16::<LittleEndian>().map(Int16::from)?;

        match value {
            Int16::Value(n) => values.push(Some(Value::from(i32::from(n)))),
            Int16::Missing => values.push(None),
            _ => todo!("unhandled i16 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int16_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        src.read_i16_into::<LittleEndian>(&mut buf)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int16::from)
            .filter_map(|value| match value {
                Int16::Value(n) => Some(Some(i32::from(n))),
                Int16::Missing => Some(None),
                Int16::EndOfVector => None,
                _ => todo!("unhandled i16 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_int32_values(src: &mut &[u8], sample_count: usize) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = src.read_i32::<LittleEndian>().map(Int32::from)?;

        match value {
            Int32::Value(n) => values.push(Some(Value::from(n))),
            Int32::Missing => values.push(None),
            _ => todo!("unhandled i32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_int32_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0; len];
        src.read_i32_into::<LittleEndian>(&mut buf)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Int32::from)
            .filter_map(|value| match value {
                Int32::Value(n) => Some(Some(n)),
                Int32::Missing => Some(None),
                Int32::EndOfVector => None,
                _ => todo!("unhandled i32 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_float_values(src: &mut &[u8], sample_count: usize) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let value = src.read_f32::<LittleEndian>().map(Float::from)?;

        match value {
            Float::Value(n) => values.push(Some(Value::from(n))),
            Float::Missing => values.push(None),
            _ => todo!("unhandled f32 value: {:?}", value),
        }
    }

    Ok(values)
}

fn read_float_array_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    for _ in 0..sample_count {
        let mut buf = vec![0.0; len];
        src.read_f32_into::<LittleEndian>(&mut buf)?;

        let vs: Vec<_> = buf
            .into_iter()
            .map(Float::from)
            .filter_map(|value| match value {
                Float::Value(n) => Some(Some(n)),
                Float::Missing => Some(None),
                Float::EndOfVector => None,
                _ => todo!("unhandled f32 array value: {:?}", value),
            })
            .collect();

        if vs.len() == 1 && vs[0].is_none() {
            values.push(None);
        } else {
            values.push(Some(Value::from(vs)));
        }
    }

    Ok(values)
}

fn read_string_values(
    src: &mut &[u8],
    sample_count: usize,
    len: usize,
) -> io::Result<Vec<Option<Value>>> {
    const NUL: u8 = 0x00;

    let mut values = Vec::with_capacity(sample_count);
    let mut buf = vec![0; len];

    for _ in 0..sample_count {
        src.read_exact(&mut buf)?;

        let data = match buf.iter().position(|&b| b == NUL) {
            Some(i) => &buf[..i],
            None => &buf[..],
        };

        let s = str::from_utf8(data).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let value = Value::from(s);

        values.push(Some(value));
    }

    Ok(values)
}

pub(super) fn read_genotype_values(
    src: &mut &[u8],
    sample_count: usize,
) -> io::Result<Vec<Option<Value>>> {
    let mut values = Vec::with_capacity(sample_count);

    match read_type(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
        Some(Type::Int8(len)) => match len {
            0 => values.push(None),
            1 => {
                for _ in 0..sample_count {
                    let value = src
                        .read_i8()
                        .map(|v| parse_genotype_values(&[v]))
                        .map(Value::from)?;

                    values.push(Some(value));
                }
            }
            _ => {
                let mut buf = vec![0; len];

                for _ in 0..sample_count {
                    src.read_i8_into(&mut buf)?;
                    let value = Value::from(parse_genotype_values(&buf));
                    values.push(Some(value));
                }
            }
        },
        ty => todo!("unhandled type: {:?}", ty),
    }

    Ok(values)
}

fn parse_genotype_values(values: &[i8]) -> String {
    use std::fmt::Write;

    let mut genotype = String::new();

    for (i, &value) in values.iter().enumerate() {
        if let Int8::EndOfVector = Int8::from(value) {
            break;
        }

        let j = (value >> 1) - 1;
        let is_phased = value & 0x01 == 1;

        if i > 0 {
            if is_phased {
                genotype.push('|');
            } else {
                genotype.push('/');
            }
        }

        if j == -1 {
            genotype.push('.');
        } else {
            let _ = write!(genotype, "{j}");
        }
    }

    genotype
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_values_with_int8_values() -> io::Result<()> {
        let mut src = &[
            0x11, // Some(Type::Int8(1))
            0x05, // Some(5)
            0x08, // Some(8)
            0x80, // None
        ][..];

        let actual = read_values(&mut src, 3)?;
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_int8_array_values() -> io::Result<()> {
        let mut src = &[
            0x21, // Some(Type::Int8(2))
            0x05, 0x08, // Some([Some(5), Some(8)])
            0x0d, 0x80, // Some([Some(13), None])
            0x15, 0x81, // Some([Some(21)])
            0x80, 0x81, // None
        ][..];

        let actual = read_values(&mut src, 4)?;
        let expected = vec![
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_int16_values() -> io::Result<()> {
        let mut src = &[
            0x12, // Some(Type::Int16(1))
            0x05, 0x00, // Some(5)
            0x08, 0x00, // Some(8)
            0x00, 0x80, // None
        ][..];

        let actual = read_values(&mut src, 3)?;
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_int16_array_values() -> io::Result<()> {
        let mut src = &[
            0x22, // Some(Type::Int16(2))
            0x05, 0x00, 0x08, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x01, 0x80, // Some([Some(21)])
            0x00, 0x80, 0x01, 0x80, // None
        ][..];

        let actual = read_values(&mut src, 4)?;
        let expected = vec![
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_int32_values() -> io::Result<()> {
        let mut src = &[
            0x13, // Some(Type::Int32(1))
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
        ][..];

        let actual = read_values(&mut src, 3)?;
        let expected = vec![Some(Value::from(5)), Some(Value::from(8)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_int32_array_values() -> io::Result<()> {
        let mut src = &[
            0x23, // Some(Type::Int32(2))
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(21)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ][..];

        let actual = read_values(&mut src, 4)?;
        let expected = vec![
            Some(Value::from(vec![Some(5), Some(8)])),
            Some(Value::from(vec![Some(13), None])),
            Some(Value::from(vec![Some(21)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_float_values() -> io::Result<()> {
        let mut src = &[
            0x15, // Some(Type::Float(1))
            0x00, 0x00, 0x00, 0x00, // Some(0.0)
            0x00, 0x00, 0x80, 0x3f, // Some(1.0)
            0x01, 0x00, 0x80, 0x7f, // None
        ][..];

        let actual = read_values(&mut src, 3)?;
        let expected = vec![Some(Value::from(0.0)), Some(Value::from(1.0)), None];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_float_array_values() -> io::Result<()> {
        let mut src = &[
            0x25, // Some(Type::Float(2))
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // Some([Some(0.0), Some(1.0)])
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // Some([Some(0.0), None])
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // Some([Some(0.0)])
            0x01, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // None
        ][..];

        let actual = read_values(&mut src, 4)?;
        let expected = vec![
            Some(Value::from(vec![Some(0.0), Some(1.0)])),
            Some(Value::from(vec![Some(0.0), None])),
            Some(Value::from(vec![Some(0.0)])),
            None,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_read_values_with_string_values() -> io::Result<()> {
        let mut src = &[
            0x47, // Some(Type::String(4))
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b'd', b'l', 0x00, // "ndl"
            b'n', b'd', b'l', b's', // "ndls"
        ][..];

        let actual = read_values(&mut src, 3)?;
        let expected = vec![
            Some(Value::from("n")),
            Some(Value::from("ndl")),
            Some(Value::from("ndls")),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_genotype_genotype_field_values() {
        // Examples from § 6.3.3 Type encoding (2021-05-13)

        assert_eq!(parse_genotype_values(&[0x02, 0x02]), "0/0");
        assert_eq!(parse_genotype_values(&[0x02, 0x04]), "0/1");
        assert_eq!(parse_genotype_values(&[0x04, 0x04]), "1/1");
        assert_eq!(parse_genotype_values(&[0x02, 0x05]), "0|1");
        assert_eq!(parse_genotype_values(&[0x00, 0x00]), "./.");
        assert_eq!(parse_genotype_values(&[0x02]), "0");
        assert_eq!(parse_genotype_values(&[0x04]), "1");
        assert_eq!(parse_genotype_values(&[0x02, 0x04, 0x06]), "0/1/2");
        assert_eq!(parse_genotype_values(&[0x02, 0x04, 0x07]), "0/1|2");
        assert_eq!(
            parse_genotype_values(&[0x02, i8::from(Int8::EndOfVector)]),
            "0"
        );
    }
}
