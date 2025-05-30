use std::{io, iter, mem};

use noodles_vcf as vcf;

use super::value::{Type, read_type};

/// BCF record filters.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filters<'r>(&'r [u8]);

impl<'r> Filters<'r> {
    pub(super) fn new(src: &'r [u8]) -> Self {
        Self(src)
    }

    fn indices(&self) -> Box<dyn Iterator<Item = io::Result<usize>> + '_> {
        fn invalid_value_error() -> io::Error {
            io::Error::new(io::ErrorKind::InvalidData, "invalid value")
        }

        let mut src = self.as_ref();

        // SAFETY: The type is guaranteed to be an integer type.
        let iter: Box<dyn Iterator<Item = io::Result<usize>>> = match read_type(&mut src).unwrap() {
            None => Box::new(iter::empty()),
            Some(Type::Int8(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i8).map_err(|_| invalid_value_error())),
            ),
            Some(Type::Int16(_)) => Box::new(src.chunks(mem::size_of::<i16>()).map(|chunk| {
                let buf = chunk.try_into().map_err(|_| invalid_value_error())?;
                usize::try_from(i16::from_le_bytes(buf)).map_err(|_| invalid_value_error())
            })),
            Some(Type::Int32(_)) => Box::new(src.chunks(mem::size_of::<i32>()).map(|chunk| {
                let buf = chunk.try_into().map_err(|_| invalid_value_error())?;
                usize::try_from(i32::from_le_bytes(buf)).map_err(|_| invalid_value_error())
            })),
            _ => unreachable!(),
        };

        iter
    }
}

impl AsRef<[u8]> for Filters<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl vcf::variant::record::Filters for Filters<'_> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn len(&self) -> usize {
        let mut src = self.as_ref();

        // SAFETY: The type is guaranteed to be an integer type.
        match read_type(&mut src).unwrap() {
            None => 0,
            Some(Type::Int8(len)) => len,
            Some(Type::Int16(len)) => len,
            Some(Type::Int32(len)) => len,
            _ => unreachable!(),
        }
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        Box::new(self.indices().map(|result| {
            result.and_then(|i| {
                header.string_maps().strings().get_index(i).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid string map index: {i}"),
                    )
                })
            })
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_indices() -> io::Result<()> {
        fn t(src: &[u8], expected: &[usize]) -> io::Result<()> {
            let filters = Filters(src);
            let actual = filters.indices().collect::<io::Result<Vec<usize>>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        // None
        t(&[0x00], &[])?;

        // Some(Type::Int8(_))
        t(&[0x11, 0x00], &[0])?;
        t(&[0x21, 0x00, 0x01], &[0, 1])?;

        // Some(Type::Int16(_))
        t(&[0x12, 0x80, 0x00], &[128])?;
        t(&[0x22, 0x00, 0x00, 0x80, 0x00], &[0, 128])?;

        // Some(Type::Int32(_))
        t(&[0x13, 0x00, 0x80, 0x00, 0x00], &[32768])?;
        t(
            &[0x23, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00],
            &[0, 32768],
        )?;

        Ok(())
    }
}
