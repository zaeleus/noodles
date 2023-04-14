use std::io::{self, Read};

use crate::lazy::record::{
    value::{Int16, Int32, Int8},
    Value,
};

use super::value::read_value;

pub fn read_string_map_index<R>(reader: &mut R) -> io::Result<usize>
where
    R: Read,
{
    let i = match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => i32::from(i),
        Some(Value::Int16(Some(Int16::Value(i)))) => i32::from(i),
        Some(Value::Int32(Some(Int32::Value(i)))) => i,
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected {{Int8, Int16, Int32}}, got {v:?}"),
            ))
        }
    };

    usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

pub fn read_string_map_indices<R>(reader: &mut R) -> io::Result<Vec<usize>>
where
    R: Read,
{
    let indices = match read_value(reader)? {
        Some(Value::Int8(Some(Int8::Value(i)))) => vec![i32::from(i)],
        Some(Value::Int8Array(indices)) => indices.into_iter().map(i32::from).collect(),
        Some(Value::Int16(Some(Int16::Value(i)))) => vec![i32::from(i)],
        Some(Value::Int16Array(indices)) => indices.into_iter().map(i32::from).collect(),
        Some(Value::Int32(Some(Int32::Value(i)))) => vec![i],
        Some(Value::Int32Array(indices)) => indices,
        None => Vec::new(),
        v => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "expected {{Int8, Int8Array, Int16, Int16Array, Int32, Int32Array}}, got {v:?}"
                ),
            ))
        }
    };

    indices
        .into_iter()
        .map(|i| usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_string_map_index() -> io::Result<()> {
        // Some(Type::Int8(Some(Int8::Value(8))))
        let data = [0x11, 0x08];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 8);

        // Some(Type::Int16(Some(Int16::Value(13))))
        let data = [0x12, 0x0d, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 13);

        // Some(Type::Int32(Some(Int32::Value(21))))
        let data = [0x13, 0x15, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 21);

        // Some(Type::String(Some(String::from("n"))))
        let data = [0x17, b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_index(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }

    #[test]
    fn test_read_string_map_indices() -> io::Result<()> {
        fn t(data: &[u8], expected: &[usize]) -> io::Result<()> {
            let mut reader = data;
            let actual = read_string_map_indices(&mut reader)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        // None
        t(&[0x00], &[])?;

        // Some(Type::Int8(Some(Int8::Value(2))))
        t(&[0x11, 0x02], &[2])?;
        // Some(Type::Int8Array(vec![2, 3]))
        t(&[0x21, 0x02, 0x03], &[2, 3])?;

        // Some(Type::Int16(Some(Int16::Value(5))))
        t(&[0x12, 0x05, 0x00], &[5])?;
        // Some(Type::Int16Array(vec![5, 8]))
        t(&[0x22, 0x05, 0x00, 0x08, 0x00], &[5, 8])?;

        // Some(Type::Int32(Some(Int32::Value(13))))
        t(&[0x13, 0x0d, 0x00, 0x00, 0x00], &[13])?;
        // Some(Type::Int32Array(vec![13, 21]))
        t(
            &[0x23, 0x0d, 0x00, 0x00, 0x00, 0x15, 0x00, 0x00, 0x00],
            &[13, 21],
        )?;

        // Some(Type::Int8(None))
        let data = [0x01];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_indices(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        // Some(Type::Int16(None))
        let data = [0x02];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_indices(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        // Some(Type::Int32(None))
        let data = [0x03];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_indices(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        // Some(Type::String(Some(String::from("n"))))
        let data = [0x17, b'n'];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_indices(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
