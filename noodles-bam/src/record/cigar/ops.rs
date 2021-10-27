use std::{io, slice};

use super::Op;

/// An iterator over the operations of a CIGAR.
///
/// This is created by calling [`super::Cigar::ops`].
pub struct Ops<'a>(slice::Iter<'a, u32>);

impl<'a> Ops<'a> {
    pub(crate) fn new(cigar: &'a [u32]) -> Self {
        Self(cigar.iter())
    }
}

impl<'a> Iterator for Ops<'a> {
    type Item = io::Result<Op>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0
            .next()
            .copied()
            .map(|n| Op::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
    }
}
