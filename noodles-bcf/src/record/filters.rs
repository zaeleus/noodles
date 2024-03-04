use std::{io, iter};

use noodles_vcf as vcf;

use super::value::{read_type, Type};

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
            Some(Type::Int16(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i16).map_err(|_| invalid_value_error())),
            ),
            Some(Type::Int32(_)) => Box::new(
                src.iter()
                    .map(|&n| usize::try_from(n as i32).map_err(|_| invalid_value_error())),
            ),
            _ => unreachable!(),
        };

        iter
    }
}

impl<'r> AsRef<[u8]> for Filters<'r> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'r> vcf::variant::record::Filters for Filters<'r> {
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
