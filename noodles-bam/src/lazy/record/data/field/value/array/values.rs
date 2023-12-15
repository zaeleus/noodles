use std::{io, marker::PhantomData, mem};

#[derive(Debug, PartialEq)]
pub struct Values<'a, N> {
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

impl<'a> Values<'a, i8> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i8>> + '_ {
        self.src.iter().map(|&b| Ok(b as i8))
    }
}

impl<'a> Values<'a, u8> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<u8>> + '_ {
        self.src.iter().copied().map(Ok)
    }
}

impl<'a> Values<'a, i16> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i16>> + '_ {
        self.src.chunks(mem::size_of::<i16>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(i16::from_le_bytes(buf))
        })
    }
}

impl<'a> Values<'a, u16> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<u16>> + '_ {
        self.src.chunks(mem::size_of::<u16>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(u16::from_le_bytes(buf))
        })
    }
}

impl<'a> Values<'a, i32> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<i32>> + '_ {
        self.src.chunks(mem::size_of::<i32>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(i32::from_le_bytes(buf))
        })
    }
}

impl<'a> Values<'a, f32> {
    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<f32>> + '_ {
        self.src.chunks(mem::size_of::<f32>()).map(|chunk| {
            let buf = chunk
                .try_into()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(f32::from_le_bytes(buf))
        })
    }
}
