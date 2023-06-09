use std::{
    cmp,
    io::{self, Write},
};

use crate::lazy::record::{
    value::{Array, Int16, Int32, Int8},
    Value,
};

use super::value::write_value;

pub fn write_string_map_index<W>(writer: &mut W, i: usize) -> io::Result<()>
where
    W: Write,
{
    if let Ok(j) = i8::try_from(i) {
        write_value(writer, Some(Value::Int8(Some(Int8::Value(j)))))
    } else if let Ok(j) = i16::try_from(i) {
        write_value(writer, Some(Value::Int16(Some(Int16::Value(j)))))
    } else if let Ok(j) = i32::try_from(i) {
        write_value(writer, Some(Value::Int32(Some(Int32::Value(j)))))
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid index: {i}"),
        ))
    }
}

pub fn write_string_map_indices<W>(writer: &mut W, indices: &[usize]) -> io::Result<()>
where
    W: Write,
{
    match indices.len() {
        0 => write_value(writer, None),
        1 => match indices.first().copied() {
            Some(i) => write_string_map_index(writer, i),
            None => unreachable!(),
        },
        _ => {
            let mut max = i32::MIN;

            for &i in indices {
                let j =
                    i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                max = cmp::max(max, j);
            }

            if max <= i32::from(Int8::MAX_VALUE) {
                let is = indices.iter().map(|&i| i as i8).collect();
                write_value(writer, Some(Value::Array(Array::Int8(is))))
            } else if max <= i32::from(Int16::MAX_VALUE) {
                let is = indices.iter().map(|&i| i as i16).collect();
                write_value(writer, Some(Value::Array(Array::Int16(is))))
            } else {
                let is = indices.iter().map(|&i| i as i32).collect();
                write_value(writer, Some(Value::Array(Array::Int32(is))))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_string_map_index() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, i: usize, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_string_map_index(buf, i)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        // Int8
        t(&mut buf, 0, &[0x11, 0x00])?;
        t(&mut buf, 127, &[0x11, 0x7f])?;

        // Int16
        t(&mut buf, 128, &[0x12, 0x80, 0x00])?;
        t(&mut buf, 32767, &[0x12, 0xff, 0x7f])?;

        // Int32
        t(&mut buf, 32768, &[0x13, 0x00, 0x80, 0x00, 0x00])?;
        t(&mut buf, 2147483647, &[0x13, 0xff, 0xff, 0xff, 0x7f])?;

        buf.clear();
        assert!(matches!(
            write_string_map_index(&mut buf, 2147483648),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_string_map_indices() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, indices: &[usize], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_string_map_indices(buf, indices)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &[], &[0x00])?;

        t(&mut buf, &[0], &[0x11, 0x00])?;
        t(&mut buf, &[128], &[0x12, 0x80, 0x00])?;
        t(&mut buf, &[32768], &[0x13, 0x00, 0x80, 0x00, 0x00])?;

        t(&mut buf, &[0, 127], &[0x21, 0x00, 0x7f])?;
        t(&mut buf, &[0, 32767], &[0x22, 0x00, 0x00, 0xff, 0x7f])?;
        t(
            &mut buf,
            &[0, 2147483647],
            &[0x23, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0x7f],
        )?;

        Ok(())
    }
}
