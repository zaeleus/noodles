use std::{
    convert::TryFrom,
    io::{self, Write},
};

use noodles_vcf as vcf;

use crate::{
    header::StringMap,
    record::value::{Float, Int16, Int32, Int8, Value},
    writer::{string_map::write_string_map_index, value::write_value},
};

pub fn write_info<W>(
    writer: &mut W,
    string_map: &StringMap,
    info: &vcf::record::Info,
) -> io::Result<()>
where
    W: Write,
{
    for field in info.values() {
        write_info_field_key(writer, string_map, field.key())?;
        write_info_field_value(writer, field.value())?;
    }

    Ok(())
}

fn write_info_field_key<W>(
    writer: &mut W,
    string_map: &StringMap,
    key: &vcf::record::info::field::Key,
) -> io::Result<()>
where
    W: Write,
{
    string_map
        .get_index_of(key.as_ref())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("info key missing from string map: {:?}", key),
            )
        })
        .and_then(|i| write_string_map_index(writer, i))
}

fn write_info_field_value<W>(
    writer: &mut W,
    value: &vcf::record::info::field::Value,
) -> io::Result<()>
where
    W: Write,
{
    use vcf::record::info::field;

    match value {
        field::Value::Integer(n) => write_info_field_integer_value(writer, *n),
        field::Value::Float(n) => write_info_field_float_value(writer, *n),
        field::Value::Flag => write_info_field_flag_value(writer),
        field::Value::Character(c) => write_info_field_character_value(writer, *c),
        field::Value::String(s) => write_info_field_string_value(writer, s),
        field::Value::IntegerArray(values) => write_info_field_integer_array_value(writer, values),
        field::Value::FloatArray(values) => write_info_field_float_array_value(writer, values),
        field::Value::CharacterArray(values) => {
            write_info_field_character_array_value(writer, values)
        }
        field::Value::StringArray(values) => write_info_field_string_array_value(writer, values),
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
            format!("invalid info field integer value: {}", n),
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
    write_value(writer, Some(Value::String(Some(c.into()))))
}

fn write_info_field_string_value<W>(writer: &mut W, s: &str) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, Some(Value::String(Some(s.into()))))
}

fn write_info_field_integer_array_value<W>(writer: &mut W, values: &[i32]) -> io::Result<()>
where
    W: Write,
{
    let max = values
        .iter()
        .max()
        .copied()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidInput))?;

    if i8::try_from(max).is_ok() {
        write_info_field_int8_array_value(writer, values)
    } else if i16::try_from(max).is_ok() {
        write_info_field_int16_array_value(writer, values)
    } else {
        write_info_field_int32_array_value(writer, values)
    }
}

fn write_info_field_int8_array_value<W>(writer: &mut W, values: &[i32]) -> io::Result<()>
where
    W: Write,
{
    let v = values
        .iter()
        .map(|&n| i8::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)))
        .collect::<Result<_, _>>()?;

    write_value(writer, Some(Value::Int8Array(v)))
}

fn write_info_field_int16_array_value<W>(writer: &mut W, values: &[i32]) -> io::Result<()>
where
    W: Write,
{
    let v = values
        .iter()
        .map(|&n| i16::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)))
        .collect::<Result<_, _>>()?;

    write_value(writer, Some(Value::Int16Array(v)))
}

fn write_info_field_int32_array_value<W>(writer: &mut W, values: &[i32]) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, Some(Value::Int32Array(values.into())))
}

fn write_info_field_float_array_value<W>(writer: &mut W, values: &[f32]) -> io::Result<()>
where
    W: Write,
{
    write_value(writer, Some(Value::FloatArray(values.into())))
}

fn write_info_field_character_array_value<W>(writer: &mut W, values: &[char]) -> io::Result<()>
where
    W: Write,
{
    let mut s = String::new();

    for (i, &c) in values.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }

        s.push(c);
    }

    write_value(writer, Some(Value::String(Some(s))))
}

fn write_info_field_string_array_value<W>(writer: &mut W, values: &[String]) -> io::Result<()>
where
    W: Write,
{
    let mut s = String::new();

    for (i, t) in values.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }

        s.push_str(t);
    }

    write_value(writer, Some(Value::String(Some(s))))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_write_info_field_value_with_integer_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: &field::Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::Integer(-2147483641);
        buf.clear();
        assert!(matches!(
            write_info_field_value(&mut buf, &value),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        let value = field::Value::Integer(-2147483640);
        t(&mut buf, &value, &[0x13, 0x08, 0x00, 0x00, 0x80])?;

        let value = field::Value::Integer(-32761);
        t(&mut buf, &value, &[0x13, 0x07, 0x80, 0xff, 0xff])?;

        let value = field::Value::Integer(-32760);
        t(&mut buf, &value, &[0x12, 0x08, 0x80])?;

        let value = field::Value::Integer(-121);
        t(&mut buf, &value, &[0x12, 0x87, 0xff])?;

        let value = field::Value::Integer(-120);
        t(&mut buf, &value, &[0x11, 0x88])?;

        let value = field::Value::Integer(0);
        t(&mut buf, &value, &[0x11, 0x00])?;

        let value = field::Value::Integer(127);
        t(&mut buf, &value, &[0x11, 0x7f])?;

        let value = field::Value::Integer(128);
        t(&mut buf, &value, &[0x12, 0x80, 0x00])?;

        let value = field::Value::Integer(32767);
        t(&mut buf, &value, &[0x12, 0xff, 0x7f])?;

        let value = field::Value::Integer(32768);
        t(&mut buf, &value, &[0x13, 0x00, 0x80, 0x00, 0x00])?;

        let value = field::Value::Integer(2147483647);
        t(&mut buf, &value, &[0x13, 0xff, 0xff, 0xff, 0x7f])?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_float_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::Float(0.0);
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x15, 0x00, 0x00, 0x00, 0x00];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_flag_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::Flag;
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x00];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_character_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::Character('n');
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x17, 0x6e];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_string_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::String(String::from("ndls"));
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x47, 0x6e, 0x64, 0x6c, 0x73];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_integer_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        fn t(buf: &mut Vec<u8>, value: &field::Value, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_info_field_value(buf, value)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        let value = field::Value::IntegerArray(vec![8, 13]);
        t(&mut buf, &value, &[0x21, 0x08, 0x0d])?;

        let value = field::Value::IntegerArray(vec![144, 233]);
        t(&mut buf, &value, &[0x22, 0x90, 0x00, 0xe9, 0x00])?;

        let value = field::Value::IntegerArray(vec![46368, 75025]);
        t(
            &mut buf,
            &value,
            &[0x23, 0x20, 0xb5, 0x00, 0x00, 0x11, 0x25, 0x01, 0x00],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_float_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::FloatArray(vec![0.0, 1.0]);
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_character_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::CharacterArray(vec!['n', 'd', 'l', 's']);
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x77, 0x6e, 0x2c, 0x64, 0x2c, 0x6c, 0x2c, 0x73];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_write_info_field_value_with_string_array_value() -> io::Result<()> {
        use vcf::record::info::field;

        let mut buf = Vec::new();
        let value = field::Value::StringArray(vec![String::from("nd"), String::from("ls")]);
        write_info_field_value(&mut buf, &value)?;

        let expected = [0x57, 0x6e, 0x64, 0x2c, 0x6c, 0x73];

        assert_eq!(buf, expected);

        Ok(())
    }
}
