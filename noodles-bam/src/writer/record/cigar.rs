use std::io;

use bytes::BufMut;
use noodles_sam::{self as sam, record::cigar::Op};

pub fn put_cigar<B>(dst: &mut B, cigar: &sam::record::Cigar) -> io::Result<()>
where
    B: BufMut,
{
    for &op in cigar.as_ref() {
        let n = encode_op(op)?;
        dst.put_u32_le(n);
    }

    Ok(())
}

fn encode_op(op: Op) -> io::Result<u32> {
    const MAX_LENGTH: u32 = (1 << 28) - 1;

    let len =
        u32::try_from(op.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    if len <= MAX_LENGTH {
        let k = op.kind() as u32;
        Ok(len << 4 | k)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CIGAR op length",
        ))
    }
}
