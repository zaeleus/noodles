use std::io;

use crate::alignment::record::cigar::Op;

/// An iterator adapter that simplifies CIGAR operations.
pub struct TrySimplify<I> {
    iter: I,
    prev_op: Option<Op>,
}

impl<I> TrySimplify<I>
where
    I: Iterator<Item = io::Result<Op>>,
{
    pub fn new(iter: I) -> Self {
        Self {
            iter,
            prev_op: None,
        }
    }
}

impl<I> Iterator for TrySimplify<I>
where
    I: Iterator<Item = io::Result<Op>>,
{
    type Item = io::Result<Op>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let Some(result) = self.iter.next() else {
                return self.prev_op.take().map(Ok);
            };

            let op = match result {
                Ok(op) => op,
                Err(e) => return Some(Err(e)),
            };

            if let Some(prev_op) = self.prev_op.replace(op) {
                if prev_op.kind() == op.kind() {
                    self.prev_op = Some(Op::new(prev_op.kind(), prev_op.len() + op.len()));
                } else {
                    return Some(Ok(prev_op));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        alignment::record::{cigar::op::Kind, Cigar},
        alignment::record_buf::Cigar as CigarBuf,
    };

    #[test]
    fn test_next() -> io::Result<()> {
        let cigar: CigarBuf = [
            Op::new(Kind::Match, 1),
            Op::new(Kind::Insertion, 2),
            Op::new(Kind::Match, 3),
            Op::new(Kind::Match, 4),
            Op::new(Kind::SoftClip, 5),
        ]
        .into_iter()
        .collect();

        let iter = TrySimplify::new(cigar.iter());
        let actual: Vec<_> = iter.collect::<Result<_, _>>()?;

        assert_eq!(
            actual,
            [
                Op::new(Kind::Match, 1),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 7),
                Op::new(Kind::SoftClip, 5)
            ]
        );

        Ok(())
    }
}
