mod op;

pub use self::op::Op;

use std::{convert::TryFrom, fmt, mem, ops::Deref};

use noodles_sam::cigar::op::Kind;

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

    pub fn mapped_len(&self) -> u32 {
        self.ops()
            .filter_map(|op| match op.kind() {
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SeqMatch | Kind::SeqMismatch => {
                    Some(op.len())
                }
                _ => None,
            })
            .sum()
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

            let data = &self.cigar[start..end];
            let op = Op::try_from(data).unwrap();

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
        let cigar = Cigar::new(&bytes);
        let mut ops = cigar.ops();
        assert_eq!(ops.next(), Some(Op::try_from(0x240).unwrap()));
        assert_eq!(ops.next(), Some(Op::try_from(0x362).unwrap()));
        assert_eq!(ops.next(), None);
    }
}
