use std::{io, mem, ops::Range};

use bytes::Buf;

use noodles_sam::{self as sam, record::data::field::value::Type};

/// A container to store the end positions of each field.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct Bounds(Vec<usize>);

impl Bounds {
    pub fn update<B>(&mut self, mut buf: B) -> io::Result<()>
    where
        B: Buf,
    {
        self.clear();

        let len = buf.remaining();

        while buf.has_remaining() {
            advance_tag(&mut buf)?;
            let ty = read_ty(&mut buf)?;
            advance_value(&mut buf, ty)?;

            let i = len - buf.remaining();
            self.push(i);
        }

        Ok(())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn get(&self, i: usize) -> Option<Range<usize>> {
        let start = i
            .checked_sub(1)
            .and_then(|j| self.0.get(j).copied())
            .unwrap_or_default();

        let end = self.0.get(i).copied()?;

        Some(start..end)
    }

    pub(super) fn push(&mut self, i: usize) {
        self.0.push(i);
    }

    pub(super) fn clear(&mut self) {
        self.0.clear();
    }
}

impl AsRef<[usize]> for Bounds {
    fn as_ref(&self) -> &[usize] {
        &self.0
    }
}

fn advance_tag<B>(buf: &mut B) -> io::Result<()>
where
    B: Buf,
{
    const LENGTH: usize = 2;

    if buf.remaining() < LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    buf.advance(LENGTH);

    Ok(())
}

fn read_ty<B>(buf: &mut B) -> io::Result<Type>
where
    B: Buf,
{
    const LENGTH: usize = mem::size_of::<u8>();

    if buf.remaining() < LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    Type::try_from(buf.get_u8()).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn advance_value<B>(buf: &mut B, ty: Type) -> io::Result<()>
where
    B: Buf,
{
    let len = match ty {
        Type::Char | Type::UInt8 => mem::size_of::<u8>(),
        Type::Int8 => mem::size_of::<i8>(),
        Type::Int16 => mem::size_of::<i16>(),
        Type::UInt16 => mem::size_of::<u16>(),
        Type::Int32 => mem::size_of::<i32>(),
        Type::UInt32 => mem::size_of::<u32>(),
        Type::Float => mem::size_of::<f32>(),
        Type::String | Type::Hex => size_of_string_value(buf)?,
        Type::Array => size_of_array_value(buf)?,
    };

    buf.advance(len);

    Ok(())
}

fn size_of_string_value<B>(buf: &mut B) -> io::Result<usize>
where
    B: Buf,
{
    const NUL: u8 = b'\0';

    buf.chunk()
        .iter()
        .position(|&b| b == NUL)
        .map(|i| i + 1)
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))
}

fn size_of_array_value<B>(buf: &mut B) -> io::Result<usize>
where
    B: Buf,
{
    use sam::record::data::field::value::Subtype;

    let subtype = Subtype::try_from(buf.get_u8())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let n = buf.get_u32_le() as usize;

    let size = match subtype {
        Subtype::Int8 => mem::size_of::<i8>(),
        Subtype::UInt8 => mem::size_of::<u8>(),
        Subtype::Int16 => mem::size_of::<i16>(),
        Subtype::UInt16 => mem::size_of::<u16>(),
        Subtype::Int32 => mem::size_of::<i32>(),
        Subtype::UInt32 => mem::size_of::<u32>(),
        Subtype::Float => mem::size_of::<f32>(),
    };

    Ok(n * size)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_update() -> io::Result<()> {
        let data = [
            b'Z', b'A', b'A', b'n', // ZA:A:n
            b'Z', b'c', b'c', 0x00, // Zc:c:0
            b'Z', b'C', b'C', 0x00, // ZC:C:0
            b'Z', b's', b's', 0x00, 0x00, // Zs:s:0
            b'Z', b'S', b'S', 0x00, 0x00, // ZS:S:0
            b'Z', b'i', b'i', 0x00, 0x00, 0x00, 0x00, // Zi:i:0
            b'Z', b'I', b'I', 0x00, 0x00, 0x00, 0x00, // ZI:I:0
            b'Z', b'f', b'f', 0x00, 0x00, 0x00, 0x00, // Zf:f:0
            b'Z', b'Z', b'Z', b'n', b'd', b'l', b's', 0x00, // ZZ:Z:ndls
            b'Z', b'H', b'H', b'C', b'A', b'F', b'E', 0x00, // ZH:H:CAFE
            b'b', b'c', b'B', b'c', 0x01, 0x00, 0x00, 0x00, 0x00, // bc:B:c,0
            b'b', b'C', b'B', b'C', 0x01, 0x00, 0x00, 0x00, 0x00, // bC:B:C,0
            b'b', b's', b'B', b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // bs:B:s,0
            b'b', b'S', b'B', b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, // bS:B:S,0
            b'b', b'i', b'B', b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, // bi:B:i,0
            b'b', b'I', b'B', b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, // bI:B:I,0
            b'b', b'f', b'B', b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, // bf:B:f,0
        ];

        let mut actual = Bounds::default();
        actual.update(&data[..])?;

        let expected = Bounds(vec![
            4, 8, 12, 17, 22, 29, 36, 43, 51, 59, 68, 77, 87, 97, 109, 121, 133,
        ]);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_get() {
        let bounds = Bounds(vec![2, 3, 5, 8]);

        assert_eq!(bounds.get(0), Some(0..2));
        assert_eq!(bounds.get(1), Some(2..3));
        assert_eq!(bounds.get(2), Some(3..5));
        assert_eq!(bounds.get(3), Some(5..8));

        assert!(bounds.get(4).is_none());
    }
}
