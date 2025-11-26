use std::{error, fmt, mem};

pub fn read_i8(src: &mut &[u8]) -> Result<i8, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(*n as i8)
}

pub fn read_i8s(src: &mut &[u8], len: usize) -> Result<Vec<i8>, DecodeError> {
    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    let values = buf.iter().map(|&n| n as i8).collect();

    *src = rest;

    Ok(values)
}

pub fn read_i16(src: &mut &[u8]) -> Result<i16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(i16::from_le_bytes(*buf))
}

pub fn read_i16s(src: &mut &[u8], len: usize) -> Result<Vec<i16>, DecodeError> {
    let len = mem::size_of::<i16>() * len;

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    let values = buf
        .chunks_exact(mem::size_of::<i16>())
        .map(|chunk| {
            // SAFETY: `chunk` is 2 bytes.
            i16::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

pub fn read_i32(src: &mut &[u8]) -> Result<i32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(i32::from_le_bytes(*buf))
}

pub fn read_i32s(src: &mut &[u8], len: usize) -> Result<Vec<i32>, DecodeError> {
    let len = mem::size_of::<i32>() * len;

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    let values = buf
        .chunks_exact(mem::size_of::<i32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            i32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

pub fn read_f32(src: &mut &[u8]) -> Result<f32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(f32::from_le_bytes(*buf))
}

pub fn read_f32s(src: &mut &[u8], len: usize) -> Result<Vec<f32>, DecodeError> {
    let len = mem::size_of::<f32>() * len;

    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    let values = buf
        .chunks_exact(mem::size_of::<f32>())
        .map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            f32::from_le_bytes(chunk.try_into().unwrap())
        })
        .collect();

    *src = rest;

    Ok(values)
}

pub fn read_string<'a>(src: &mut &'a [u8], len: usize) -> Result<&'a [u8], DecodeError> {
    let (buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    Ok(buf)
}

#[derive(Debug, Eq, PartialEq)]
pub enum DecodeError {
    UnexpectedEof,
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_i8() {
        let mut src = &[0x00][..];
        assert_eq!(read_i8(&mut src), Ok(0));

        let mut src = &[][..];
        assert_eq!(read_i8(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_i8s() {
        let mut src = &[0x00][..];
        assert_eq!(read_i8s(&mut src, 1), Ok(vec![0]));

        let mut src = &[][..];
        assert_eq!(read_i8s(&mut src, 1), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_i16() {
        let mut src = &[0x00, 0x00][..];
        assert_eq!(read_i16(&mut src), Ok(0));

        let mut src = &[][..];
        assert_eq!(read_i16(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_i16s() {
        let mut src = &[0x00, 0x00][..];
        assert_eq!(read_i16s(&mut src, 1), Ok(vec![0]));

        let mut src = &[][..];
        assert_eq!(read_i16s(&mut src, 1), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_i32() {
        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_i32(&mut src), Ok(0));

        let mut src = &[][..];
        assert_eq!(read_i32(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_i32s() {
        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_i32s(&mut src, 1), Ok(vec![0]));

        let mut src = &[][..];
        assert_eq!(read_i32s(&mut src, 1), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_f32() {
        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_f32(&mut src), Ok(0.0));

        let mut src = &[][..];
        assert_eq!(read_f32(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_f32s() {
        let mut src = &[0x00, 0x00, 0x00, 0x00][..];
        assert_eq!(read_f32s(&mut src, 1), Ok(vec![0.0]));

        let mut src = &[][..];
        assert_eq!(read_f32s(&mut src, 1), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_string() {
        let mut src = &[b'n', b'd', b'l', b's'][..];
        assert_eq!(read_string(&mut src, 4), Ok(&b"ndls"[..]));

        let mut src = &[][..];
        assert_eq!(read_string(&mut src, 4), Err(DecodeError::UnexpectedEof));
    }
}
