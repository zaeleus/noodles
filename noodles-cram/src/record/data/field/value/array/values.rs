use std::{io, marker::PhantomData, mem};

use noodles_sam as sam;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Values<'c, N> {
    src: &'c [u8],
    len: usize,
    _marker: PhantomData<N>,
}

impl<'c, N> Values<'c, N> {
    pub(crate) fn new(src: &'c [u8], len: usize) -> Self {
        Self {
            src,
            len,
            _marker: PhantomData,
        }
    }
}

const OFFSET: usize = 5;

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i8> for Values<'c, i8> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i8>> + '_> {
        Box::new(self.src[OFFSET..].iter().copied().map(|n| n as i8).map(Ok))
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u8> for Values<'c, u8> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u8>> + '_> {
        Box::new(self.src[OFFSET..].iter().copied().map(Ok))
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i16> for Values<'c, i16> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i16>> + '_> {
        Box::new(
            self.src[OFFSET..]
                .chunks_exact(mem::size_of::<i16>())
                .map(|buf| i16::from_le_bytes(buf.try_into().unwrap()))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u16> for Values<'c, u16> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u16>> + '_> {
        Box::new(
            self.src[OFFSET..]
                .chunks_exact(mem::size_of::<u16>())
                .map(|buf| u16::from_le_bytes(buf.try_into().unwrap()))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, i32> for Values<'c, i32> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i32>> + '_> {
        Box::new(
            self.src[OFFSET..]
                .chunks_exact(mem::size_of::<i32>())
                .map(|buf| i32::from_le_bytes(buf.try_into().unwrap()))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, u32> for Values<'c, u32> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<u32>> + '_> {
        Box::new(
            self.src[OFFSET..]
                .chunks_exact(mem::size_of::<u32>())
                .map(|buf| u32::from_le_bytes(buf.try_into().unwrap()))
                .map(Ok),
        )
    }
}

impl<'c> sam::alignment::record::data::field::value::array::Values<'c, f32> for Values<'c, f32> {
    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<f32>> + '_> {
        Box::new(
            self.src[OFFSET..]
                .chunks_exact(mem::size_of::<f32>())
                .map(|buf| f32::from_le_bytes(buf.try_into().unwrap()))
                .map(Ok),
        )
    }
}
