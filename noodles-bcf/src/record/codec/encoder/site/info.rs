use std::{
    cmp,
    io::{self, Write},
};

use noodles_vcf as vcf;

use crate::{
    header::string_maps::StringStringMap,
    lazy::record::{
        value::{Array, Float, Int16, Int32, Int8},
        Value,
    },
    record::codec::encoder::{string_map::write_string_map_index, value::write_value},
};

const MISSING_VALUE: char = '.';
const DELIMITER: char = ',';

pub fn write_info<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    info: &vcf::record::Info,
) -> io::Result<()>
where
    W: Write,
{
    for (key, value) in info.as_ref() {
        write_info_field(writer, string_string_map, key, value.as_ref())?;
    }

    Ok(())
}

fn write_info_field<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    key: &vcf::record::info::field::Key,
    value: Option<&vcf::record::info::field::Value>,
) -> io::Result<()>
where
    W: Write,
{
    write_info_field_key(writer, string_string_map, key)?;
    write_info_field_value(writer, value)?;
    Ok(())
}

fn write_info_field_key<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    key: &vcf::record::info::field::Key,
) -> io::Result<()>
where
    W: Write,
{
    string_string_map
        .get_index_of(key.as_ref())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("info key missing from string map: {key:?}"),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}

fn write_info_field_value<W>(
    writer: &mut W,
    value: Option<&vcf::record::info::field::Value>,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::info::field;

    match value {
        Some(field::Value::Integer(n)) => write_info_field_integer_value(writer, *n),
        Some(field::Value::Float(n)) => write_info_field_float_value(writer, *n),
        Some(field::Value::Flag) => write_info_field_flag_value(writer),
        Some(field::Value::Character(c)) => write_info_field_character_value(writer, *c),
        Some(field::Value::String(s)) => write_info_field_string_value(writer, s),
        Some(field::Value::Array(field::value::Array::Integer(values))) => {
            write_info_field_integer_array_value(writer, values)
        }
        Some(field::Value::Array(field::value::Array::Float(values))) => {
            write_info_field_float_array_value(writer, values)
        }
        Some(field::Value::Array(field::value::Array::Character(values))) => {
            write_info_field_character_array_value(writer, values)
        }
        Some(field::Value::Array(field::value::Array::String(values))) => {
            write_info_field_string_array_value(writer, values)
        }
        _ => todo!("unhandled INFO field value: {:?}", value),
    }
}

fn write_info_field_integer_value<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    if n >= 0 {
        if n <= i32::from(Int8::MAX_VALUE) {
            write_value(writer, Some(Value::Int8(Some(Int8::Value(n as i8)))))
        } else if n <= i32::from(Int16::MAX_VALUE) {
            write_value(writer, Some(Value::Int16(Some(Int16::Value(n as i16)))))
        } else {
            write_value(writer, Some(Value::Int32(Some(Int32::Value(n)))))
        }
    } else if n >= i32::from(Int8::MIN_VALUE) {
        write_value(writer, Some(Value::Int8(Some(Int8::Value(n as i8)))))
    } else if n >= i32::from(Int16::MIN_VALUE) {
        write_value(writer, Some(Value::Int16(Some(Int16::Value(n as i16)))))
    } else if n >= Int32::MIN_VALUE {
        write_value(writer, Some(Value::Int32(Some(Int32::Value(n)))))
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid info field integer value: {n}"),
        ))
    }
}

fn write_info_field_float_value<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, Some(Value::Float(Some(Float::Value(n)))))
}

fn write_info_field_flag_value<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, None)
}

fn write_info_field_character_value<W>(writer: &mut W, c: char) -> io::Result<()>
where
    W: Write,
{
    let s = c.to_string();
    write_value(writer, Some(Value::String(Some(&s))))
}

fn write_info_field_string_value<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, Some(Value::String(Some(s))))
}

fn write_info_field_integer_array_value<W>(writer: &mut W, values: &[Option<i32>]) -> io::Result<()>
where
    W: Write,
{
    if values.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "info field integer array cannot be empty",
        ));
    }

    let (mut min, mut max) = (i32::MAX, i32::MIN);

    for value in values {
        let n = value.unwrap_or_default();
        min = cmp::min(min, n);
        max = cmp::max(max, n);
    }

    if min >= i32::from(Int8::MIN_VALUE) {
        if max <= i32::from(Int8::MAX_VALUE) {
            write_info_field_int8_array_value(writer, values)
        } else if max <= i32::from(Int16::MAX_VALUE) {
            write_info_field_int16_array_value(writer, values)
        } else {
            write_info_field_int32_array_value(writer, values)
        }
    } else if min >= i32::from(Int16::MIN_VALUE) {
        if max <= i32::from(Int16::MAX_VALUE) {
            write_info_field_int16_array_value(writer, values)
        } else {
            write_info_field_int32_array_value(writer, values)
        }
    } else if min >= Int32::MIN_VALUE {
        write_info_field_int32_array_value(writer, values)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid info field integer array value: {min}"),
        ))
    }
}

fn write_info_field_int8_array_value<W>(writer: &mut W, values: &[Option<i32>]) -> io::Result<()>
where
    W: Write,
{
    let mut vs = Vec::with_capacity(values.len());

    for value in values {
        let v = match value {
            Some(n) => i8::try_from(*n)
                .map(Int8::from)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
            None => Int8::Missing,
        };

        let n = match v {
            Int8::Value(n) => n,
            Int8::Missing => i8::from(v),
            _ => todo!("unhandled i16 array value: {:?}", v),
        };

        vs.push(n);
    }

    write_value(writer, Some(Value::Array(Array::Int8(vs))))
}

fn write_info_field_int16_array_value<W>(writer: &mut W, values: &[Option<i32>]) -> io::Result<()>
where
    W: Write,
{
    let mut vs = Vec::with_capacity(values.len());

    for value in values {
        let v = match value {
            Some(n) => i16::try_from(*n)
                .map(Int16::from)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
            None => Int16::Missing,
        };

        let n = match v {
            Int16::Value(n) => n,
            Int16::Missing => i16::from(v),
            _ => todo!("unhandled i16 array value: {:?}", v),
        };

        vs.push(n);
    }

    write_value(writer, Some(Value::Array(Array::Int16(vs))))
}

fn write_info_field_int32_array_value<W>(writer: &mut W, values: &[Option<i32>]) -> io::Result<()>
where
    W: Write,
{
    let vs = values
        .iter()
        .map(|value| value.map(Int32::from).unwrap_or(Int32::Missing))
        .map(|value| match value {
            Int32::Value(n) => n,
            Int32::Missing => i32::from(value),
            _ => todo!("unhandled i32 array value: {:?}", value),
        })
        .collect();

    write_value(writer, Some(Value::Array(Array::Int32(vs))))
}

fn write_info_field_float_array_value<W>(writer: &mut W, values: &[Option<f32>]) -> io::Result<()>
where
    W: Write,
{
    let vs = values
        .iter()
        .map(|value| value.map(Float::from).unwrap_or(Float::Missing))
        .map(|value| match value {
            Float::Value(n) => n,
            Float::Missing => f32::from(value),
            _ => todo!("unhandled f32 array value: {:?}", value),
        })
        .collect();

    write_value(writer, Some(Value::Array(Array::Float(vs))))
}

fn write_info_field_character_array_value<W>(
    writer: &mut W,
    values: &[Option<char>],
) -> io::Result<()>
where
    W: Write,
{
    let mut s = String::new();

    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            s.push(DELIMITER);
        }

        let c = value.unwrap_or(MISSING_VALUE);
        s.push(c);
    }

    write_value(writer, Some(Value::String(Some(&s))))
}

fn write_info_field_string_array_value<W>(
    writer: &mut W,
    values: &[Option<String>],
) -> io::Result<()>
where
    W: Write,
{
    let mut s = String::new();

    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            s.push(DELIMITER);
        }

        if let Some(t) = value {
            s.push_str(t);
        } else {
            s.push(MISSING_VALUE);
        }
    }

    write_value(writer, Some(Value::String(Some(&s))))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_write_info_field_value_with_integer_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: &field::Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, Some(value))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::from(-2147483641);
        buf.clear();
        assert!(matches!(
            write_info_field_value(&mut buf, Some(&value)),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let value = field::Value::from(-2147483640);
        t(&mut buf, &value, &[0x13, 0x08, 0x00, 0x00, 0x80])?;

        let value = field::Value::from(-32761);
        t(&mut buf, &value, &[0x13, 0x07, 0x80, 0xff, 0xff])?;

        let value = field::Value::from(-32760);
        t(&mut buf, &value, &[0x12, 0x08, 0x80])?;

        let value = field::Value::from(-121);
        t(&mut buf, &value, &[0x12, 0x87, 0xff])?;

        let value = field::Value::from(-120);
        t(&mut buf, &value, &[0x11, 0x88])?;

        let value = field::Value::from(0);
        t(&mut buf, &value, &[0x11, 0x00])?;

        let value = field::Value::from(127);
        t(&mut buf, &value, &[0x11, 0x7f])?;

        let value = field::Value::from(128);
        t(&mut buf, &value, &[0x12, 0x80, 0x00])?;

        let value = field::Value::from(32767);
        t(&mut buf, &value, &[0x12, 0xff, 0x7f])?;

        let value = field::Value::from(32768);
        t(&mut buf, &value, &[0x13, 0x00, 0x80, 0x00, 0x00])?;

        let value = field::Value::from(2147483647);
        t(&mut buf, &value, &[0x13, 0xff, 0xff, 0xff, 0x7f])?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_float_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::from(0.0);
        write_info_field_value(&mut buf, Some(&value))?;

        let expected = [0x15, 0x00, 0x00, 0x00, 0x00];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_flag_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::Flag;
        write_info_field_value(&mut buf, Some(&value))?;

        let expected = [0x00];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_character_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::from('n');
        write_info_field_value(&mut buf, Some(&value))?;

        let expected = [0x17, 0x6e];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_string_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::from("ndls");
        write_info_field_value(&mut buf, Some(&value))?;

        let expected = [0x47, 0x6e, 0x64, 0x6c, 0x73];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_integer_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: Option<&field::Value>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::from(vec![Some(-2147483641), Some(-2147483640)]);
        buf.clear();
        assert!(matches!(
            write_info_field_value(&mut buf, Some(&value)),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let value = field::Value::from(vec![Some(-2147483640), Some(-2147483639)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x08, 0x00, 0x00, 0x80, 0x09, 0x00, 0x00, 0x80],
        )?;
        let value = field::Value::from(vec![Some(-2147483640), None]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x08, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00, 0x80],
        )?;

        let value = field::Value::from(vec![Some(-32761), Some(-32760)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x07, 0x80, 0xff, 0xff, 0x08, 0x80, 0xff, 0xff],
        )?;
        let value = field::Value::from(vec![Some(-32761), None]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x07, 0x80, 0xff, 0xff, 0x00, 0x00, 0x00, 0x80],
        )?;

        let value = field::Value::from(vec![Some(-32760), Some(-32759)]);
        t(&mut buf, Some(&value), &[0x22, 0x08, 0x80, 0x09, 0x80])?;
        let value = field::Value::from(vec![Some(-32760), None]);
        t(&mut buf, Some(&value), &[0x22, 0x08, 0x80, 0x00, 0x80])?;

        let value = field::Value::from(vec![Some(-121), Some(-120)]);
        t(&mut buf, Some(&value), &[0x22, 0x87, 0xff, 0x88, 0xff])?;
        let value = field::Value::from(vec![Some(-121), None]);
        t(&mut buf, Some(&value), &[0x22, 0x87, 0xff, 0x00, 0x80])?;

        let value = field::Value::from(vec![Some(-120), Some(-119)]);
        t(&mut buf, Some(&value), &[0x21, 0x88, 0x89])?;
        let value = field::Value::from(vec![Some(-120), None]);
        t(&mut buf, Some(&value), &[0x21, 0x88, 0x80])?;

        let value = field::Value::from(vec![None, Some(0), Some(1)]);
        t(&mut buf, Some(&value), &[0x31, 0x80, 0x00, 0x01])?;
        let value = field::Value::from(vec![Some(-1), Some(0), Some(1)]);
        t(&mut buf, Some(&value), &[0x31, 0xff, 0x00, 0x01])?;
        let value = field::Value::from(vec![Some(-1), Some(0), None]);
        t(&mut buf, Some(&value), &[0x31, 0xff, 0x00, 0x80])?;

        let value = field::Value::from(vec![Some(126), Some(127)]);
        t(&mut buf, Some(&value), &[0x21, 0x7e, 0x7f])?;
        let value = field::Value::from(vec![None, Some(127)]);
        t(&mut buf, Some(&value), &[0x21, 0x80, 0x7f])?;

        let value = field::Value::from(vec![Some(127), Some(128)]);
        t(&mut buf, Some(&value), &[0x22, 0x7f, 0x00, 0x80, 0x00])?;
        let value = field::Value::from(vec![None, Some(128)]);
        t(&mut buf, Some(&value), &[0x22, 0x00, 0x80, 0x80, 0x00])?;

        let value = field::Value::from(vec![Some(32766), Some(32767)]);
        t(&mut buf, Some(&value), &[0x22, 0xfe, 0x7f, 0xff, 0x7f])?;
        let value = field::Value::from(vec![None, Some(32767)]);
        t(&mut buf, Some(&value), &[0x22, 0x00, 0x80, 0xff, 0x7f])?;

        let value = field::Value::from(vec![Some(32767), Some(32768)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0xff, 0x7f, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00],
        )?;
        let value = field::Value::from(vec![None, Some(32768)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x00, 0x00, 0x00, 0x80, 0x00, 0x80, 0x00, 0x00],
        )?;

        let value = field::Value::from(vec![Some(2147483646), Some(2147483647)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0xfe, 0xff, 0xff, 0x7f, 0xff, 0xff, 0xff, 0x7f],
        )?;
        let value = field::Value::from(vec![None, Some(2147483647)]);
        t(
            &mut buf,
            Some(&value),
            &[0x23, 0x00, 0x00, 0x00, 0x80, 0xff, 0xff, 0xff, 0x7f],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_float_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: Option<&field::Value>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::from(vec![Some(0.0), Some(1.0)]);
        t(
            &mut buf,
            Some(&value),
            &[0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f],
        )?;

        let value = field::Value::from(vec![Some(0.0), None]);
        t(
            &mut buf,
            Some(&value),
            &[0x25, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_character_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: Option<&field::Value>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::from(vec![Some('n'), Some('d'), Some('l'), Some('s')]);
        t(
            &mut buf,
            Some(&value),
            &[0x77, 0x6e, 0x2c, 0x64, 0x2c, 0x6c, 0x2c, 0x73],
        )?;

        let value = field::Value::from(vec![Some('n'), Some('d'), Some('l'), None]);
        t(
            &mut buf,
            Some(&value),
            &[0x77, 0x6e, 0x2c, 0x64, 0x2c, 0x6c, 0x2c, 0x2e],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_string_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: Option<&field::Value>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::from(vec![Some(String::from("nd")), Some(String::from("ls"))]);
        t(
            &mut buf,
            Some(&value),
            &[0x57, 0x6e, 0x64, 0x2c, 0x6c, 0x73],
        )?;

        let value = field::Value::from(vec![Some(String::from("nd")), None]);
        t(&mut buf, Some(&value), &[0x47, 0x6e, 0x64, 0x2c, 0x2e])?;

        Ok(())
    }
}
