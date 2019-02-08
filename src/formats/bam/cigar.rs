pub use self::op::Op;

pub mod op;

use std::{fmt, mem, ops::Deref};

#[derive(Debug)]
pub struct Cigar<'a>(&'a [u8]);

impl<'a> Cigar<'a> {
    pub fn new(bytes: &[u8]) -> Cigar {
        Cigar(bytes)
    }

    pub fn ops(&self) -> Ops {
        Ops {
            cigar: self.0,
            i: 0,
        }
    }
}

impl<'a> fmt::Display for Cigar<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for op in self.ops() {
            write!(f, "{}", op)?;
        }

        Ok(())
    }
}

impl<'a> Deref for Cigar<'a> {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        self.0
    }
}

pub struct Ops<'a> {
    cigar: &'a [u8],
    i: usize,
}

impl<'a> Iterator for Ops<'a> {
    type Item = Op;

    fn next(&mut self) -> Option<Self::Item> {
        let size = mem::size_of::<u32>();
        let start = self.i * size;

        if start < self.cigar.len() {
            let end = start + size;
            let op = Op::from_bytes(&self.cigar[start..end]);
            self.i += 1;
            Some(op)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bytes() {
        let bytes = [0x40, 0x02, 0x00, 0x00, 0x62, 0x03, 0x00, 0x00];
        let cigar =  Cigar::new(&bytes);
        let mut ops = cigar.ops();
        assert_eq!(ops.next(), Some(Op::from_u32(0x240)));
        assert_eq!(ops.next(), Some(Op::from_u32(0x362)));
        assert_eq!(ops.next(), None);
    }
}
