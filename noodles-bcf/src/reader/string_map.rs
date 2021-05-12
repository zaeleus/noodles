use std::{
    convert::TryFrom,
    io::{self, Read},
};

use crate::record::{
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
                format!("expected {{Int8, Int16, Int32}}, got {:?}", v),
            ))
        }
    };

    usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_string_map_index() -> io::Result<()> {
        // [Type::Int8(1), 8]
        let data = [0x11, 0x08];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 8);

        // [Type::Int16(1), 13]
        let data = [0x12, 0x0d, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 13);

        // [Type::Int32(1), 21]
        let data = [0x13, 0x15, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_string_map_index(&mut reader)?, 21);

        // [Type::String(1), "n"]
        let data = [0x17, 0x6e];
        let mut reader = &data[..];
        assert!(matches!(
            read_string_map_index(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
