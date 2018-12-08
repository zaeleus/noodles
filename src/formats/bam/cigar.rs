pub use self::op::Op;

pub mod op;

use std::fmt;
use std::mem;
use std::ops::Deref;

use byteorder::{ByteOrder, LittleEndian};

#[derive(Debug)]
pub struct Cigar {
    cigar: Vec<u32>,
}

impl Cigar {
    pub fn from_bytes(bytes: &[u8]) -> Cigar {
        let cigar = bytes.chunks(mem::size_of::<u32>())
            .map(|c| LittleEndian::read_u32(c))
            .collect();

        Cigar::new(cigar)
    }

    pub fn new(cigar: Vec<u32>) -> Cigar {
        Cigar { cigar }
    }

    pub fn ops<'a>(&'a self) -> Ops<impl Iterator<Item = &'a u32>> {
        Ops(self.cigar.iter())
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for op in self.ops() {
            write!(f, "{}", op)?;
        }

        Ok(())
    }
}

impl Deref for Cigar {
    type Target = [u32];

    fn deref(&self) -> &[u32] {
        &self.cigar
    }
}

pub struct Ops<I>(I);

impl<'a, I: Iterator<Item = &'a u32>> Iterator for Ops<I> {
    type Item = Op;

    fn next(&mut self) -> Option<Op> {
        self.0.next().map(|&u| Op::from_u32(u))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_bytes() {
        let bytes = [0x40, 0x02, 0x00, 0x00, 0x62, 0x03, 0x00, 0x00];
        let cigar =  Cigar::from_bytes(&bytes);
        assert_eq!(cigar.cigar.len(), 2);
        assert_eq!(&cigar.cigar[0], &0x240);
        assert_eq!(&cigar.cigar[1], &0x362);
    }
}
