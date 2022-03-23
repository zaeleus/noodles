use std::{io, mem};

use bytes::Buf;
use noodles_sam::{
    self as sam,
    record::cigar::{op::Kind, Op},
};

pub(super) fn get_cigar<B>(
    buf: &mut B,
    cigar: &mut sam::record::Cigar,
    n_cigar_op: usize,
) -> io::Result<()>
where
    B: Buf,
{
    if buf.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_cigar_op(buf.get_u32_le())?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

fn decode_cigar_op(n: u32) -> io::Result<Op> {
    let kind = match n & 0x0f {
        0 => Kind::Match,
        1 => Kind::Insertion,
        2 => Kind::Deletion,
        3 => Kind::Skip,
        4 => Kind::SoftClip,
        5 => Kind::HardClip,
        6 => Kind::Pad,
        7 => Kind::SequenceMatch,
        8 => Kind::SequenceMismatch,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid CIGAR op kind",
            ))
        }
    };

    let len = usize::try_from(n >> 4).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(Op::new(kind, len))
}
