#![doc(hidden)]

use std::io;

use crate::alignment::record::cigar::op::Kind;

pub struct TrySimplify<I> {
    iter: I,
    prev_op: Option<(Kind, usize)>,
}

impl<I> TrySimplify<I>
where
    I: Iterator<Item = io::Result<(Kind, usize)>>,
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
    I: Iterator<Item = io::Result<(Kind, usize)>>,
{
    type Item = io::Result<(Kind, usize)>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let Some(result) = self.iter.next() else {
                return self.prev_op.take().map(Ok);
            };

            let (kind, len) = match result {
                Ok(op) => op,
                Err(e) => return Some(Err(e)),
            };

            if let Some((prev_kind, prev_len)) = self.prev_op.replace((kind, len)) {
                if prev_kind == kind {
                    self.prev_op = Some((prev_kind, prev_len + len));
                } else {
                    return Some(Ok((prev_kind, prev_len)));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        alignment::record::cigar::{op::Kind, Op},
        record::Cigar as CigarBuf,
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

        let iter = TrySimplify::new(cigar.as_ref().iter().map(|op| Ok((op.kind(), op.len()))));
        let actual: Vec<_> = iter.collect::<Result<_, _>>()?;

        assert_eq!(
            actual,
            [
                (Kind::Match, 1),
                (Kind::Insertion, 2),
                (Kind::Match, 7),
                (Kind::SoftClip, 5)
            ]
        );

        Ok(())
    }
}
