use std::{convert::TryFrom, io, mem};

use super::Op;

/// An iterator over the operations of a CIGAR.
///
/// This is created by calling [`super::Cigar::ops`].
pub struct Ops<'a> {
    cigar: &'a [u8],
    i: usize,
}

impl<'a> Ops<'a> {
    pub(crate) fn new(cigar: &'a [u8]) -> Self {
        Self { cigar, i: 0 }
    }
}

impl<'a> Iterator for Ops<'a> {
    type Item = io::Result<Op>;

    fn next(&mut self) -> Option<Self::Item> {
        let size = mem::size_of::<u32>();
        let start = self.i * size;

        if start < self.cigar.len() {
            let data = &self.cigar[start..];

            let option = Op::try_from(data)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                .transpose();

            self.i += 1;

            option
        } else {
            None
        }
    }
}
