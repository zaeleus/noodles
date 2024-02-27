use std::{io, iter};

use noodles_vcf as vcf;

use super::value::{read_type, Type};

/// BCF record filters.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Filters<'a>(&'a [u8]);

impl<'a> Filters<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any filters.
    pub fn is_empty(&self) -> io::Result<bool> {
        self.len().map(|len| len == 0)
    }

    /// Returns the number of filters.
    pub fn len(&self) -> io::Result<usize> {
        let mut src = self.as_ref();

        match read_type(&mut src)? {
            None => Ok(0),
            Some(Type::Int8(len)) => Ok(len),
            Some(Type::Int16(len)) => Ok(len),
            Some(Type::Int32(len)) => Ok(len),
            _ => Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
        }
    }

    /// Returns an iterator over filters.
    pub fn iter<'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> io::Result<Box<dyn Iterator<Item = io::Result<&'h str>> + 'a>> {
        Ok(Box::new(self.indices()?.map(|result| {
            result.and_then(|i| {
                header.string_maps().strings().get_index(i).ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid string map index: {i}"),
                    )
                })
            })
        })))
    }

    fn indices(&self) -> io::Result<Box<dyn Iterator<Item = io::Result<usize>> + '_>> {
        fn invalid_value_error() -> io::Error {
            io::Error::new(io::ErrorKind::InvalidData, "invalid value")
        }

        let mut src = self.as_ref();

        let iter: Box<dyn Iterator<Item = io::Result<usize>>> = match read_type(&mut src)? {
            None => Box::new(iter::empty()),
            Some(Type::Int8(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i8).map_err(|_| invalid_value_error())),
            ),
            Some(Type::Int16(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i16).map_err(|_| invalid_value_error())),
            ),
            Some(Type::Int32(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i32).map_err(|_| invalid_value_error())),
            ),
            _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
        };

        Ok(iter)
    }
}

impl<'a> AsRef<[u8]> for Filters<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}
