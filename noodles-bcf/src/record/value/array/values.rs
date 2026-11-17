use std::{io, marker::PhantomData};

use noodles_vcf as vcf;

use crate::record::codec::value::{Float, Int8, Int16, Int32};

pub(crate) struct Values<'a, N> {
    src: &'a [u8],
    _marker: PhantomData<N>,
}

impl<'a, N> Values<'a, N> {
    pub(crate) fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            _marker: PhantomData,
        }
    }
}

impl<N> AsRef<[u8]> for Values<'_, N> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl Values<'_, i8> {
    fn count(&self) -> usize {
        self.iter()
            .take_while(|value| !matches!(value, Int8::EndOfVector))
            .count()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Int8> + '_ {
        self.src.iter().map(|&n| Int8::from(n as i8))
    }
}

impl<'a> vcf::variant::record::info::field::value::array::Values<'a, i32> for Values<'a, i8> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int8::Value(n) => Some(Ok(Some(i32::from(n)))),
            Int8::Missing => Some(Ok(None)),
            Int8::EndOfVector => None,
            Int8::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl<'a> vcf::variant::record::samples::series::value::array::Values<'a, i32> for Values<'a, i8> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int8::Value(n) => Some(Ok(Some(i32::from(n)))),
            Int8::Missing => Some(Ok(None)),
            Int8::EndOfVector => None,
            Int8::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl Values<'_, i16> {
    fn count(&self) -> usize {
        self.iter()
            .take_while(|value| !matches!(value, Int16::EndOfVector))
            .count()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Int16> + '_ {
        let (chunks, []) = self.src.as_chunks() else {
            panic!();
        };

        chunks.iter().map(|chunk| {
            let n = i16::from_le_bytes(*chunk);
            Int16::from(n)
        })
    }
}

impl<'a> vcf::variant::record::info::field::value::array::Values<'a, i32> for Values<'a, i16> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int16::Value(n) => Some(Ok(Some(i32::from(n)))),
            Int16::Missing => Some(Ok(None)),
            Int16::EndOfVector => None,
            Int16::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl<'a> vcf::variant::record::samples::series::value::array::Values<'a, i32> for Values<'a, i16> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int16::Value(n) => Some(Ok(Some(i32::from(n)))),
            Int16::Missing => Some(Ok(None)),
            Int16::EndOfVector => None,
            Int16::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl Values<'_, i32> {
    fn count(&self) -> usize {
        self.iter()
            .take_while(|value| !matches!(value, Int32::EndOfVector))
            .count()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Int32> + '_ {
        let (chunks, []) = self.src.as_chunks() else {
            panic!();
        };

        chunks.iter().map(|chunk| {
            let n = i32::from_le_bytes(*chunk);
            Int32::from(n)
        })
    }
}

impl<'a> vcf::variant::record::info::field::value::array::Values<'a, i32> for Values<'a, i32> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int32::Value(n) => Some(Ok(Some(n))),
            Int32::Missing => Some(Ok(None)),
            Int32::EndOfVector => None,
            Int32::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl<'a> vcf::variant::record::samples::series::value::array::Values<'a, i32> for Values<'a, i32> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<i32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Int32::Value(n) => Some(Ok(Some(n))),
            Int32::Missing => Some(Ok(None)),
            Int32::EndOfVector => None,
            Int32::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl Values<'_, f32> {
    fn count(&self) -> usize {
        self.iter()
            .take_while(|value| !matches!(value, Float::EndOfVector))
            .count()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = Float> + '_ {
        let (chunks, []) = self.src.as_chunks() else {
            panic!();
        };

        chunks.iter().map(|chunk| {
            let n = f32::from_le_bytes(*chunk);
            Float::from(n)
        })
    }
}

impl<'a> vcf::variant::record::info::field::value::array::Values<'a, f32> for Values<'a, f32> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<f32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Float::Value(n) => Some(Ok(Some(n))),
            Float::Missing => Some(Ok(None)),
            Float::EndOfVector => None,
            Float::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}

impl<'a> vcf::variant::record::samples::series::value::array::Values<'a, f32> for Values<'a, f32> {
    fn len(&self) -> usize {
        self.count()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<f32>>> + '_> {
        Box::new(self.iter().map_while(|value| match value {
            Float::Value(n) => Some(Ok(Some(n))),
            Float::Missing => Some(Ok(None)),
            Float::EndOfVector => None,
            Float::Reserved(_) => Some(Err(io::Error::from(io::ErrorKind::InvalidData))),
        }))
    }
}
