use std::{marker::PhantomData, mem};

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

impl<'a> Values<'a, i8> {
    pub(crate) fn iter(&self) -> impl Iterator<Item = i8> + '_ {
        self.src.iter().map(|&n| n as i8)
    }
}

impl<'a> Values<'a, i16> {
    pub(crate) fn iter(&self) -> impl Iterator<Item = i16> + '_ {
        self.src.chunks_exact(mem::size_of::<i16>()).map(|chunk| {
            // SAFETY: `chunk` is 2 bytes.
            i16::from_le_bytes(chunk.try_into().unwrap())
        })
    }
}

impl<'a> Values<'a, i32> {
    pub(crate) fn iter(&self) -> impl Iterator<Item = i32> + '_ {
        self.src.chunks_exact(mem::size_of::<i32>()).map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            i32::from_le_bytes(chunk.try_into().unwrap())
        })
    }
}

impl<'a> Values<'a, f32> {
    pub(crate) fn iter(&self) -> impl Iterator<Item = f32> + '_ {
        self.src.chunks_exact(mem::size_of::<f32>()).map(|chunk| {
            // SAFETY: `chunk` is 4 bytes.
            f32::from_le_bytes(chunk.try_into().unwrap())
        })
    }
}
