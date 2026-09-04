mod op;

use std::{error, fmt, io};

use noodles_sam::{
    self as sam,
    alignment::record::{
        Cigar, CigarRef,
        cigar::{Op, op::Kind},
    },
};

use self::op::encode_op;
use super::num::{write_u16_le, write_u32_le};

#[derive(Debug)]
pub enum EncodeError {
    Io(io::Error),
    InvalidOp(op::EncodeError),
}

impl error::Error for EncodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidOp(e) => Some(e),
        }
    }
}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

// § 4.2.2 "N_CIGAR_OP field" (2023-11-16): "For an alignment with more [than 65535] CIGAR
// operations, BAM [...] sets `CIGAR` to `kSmN` as a placeholder, where `k` equals `l_seq`, `m` is
// the reference sequence length in the alignment..."
pub(super) fn overflowing_write_cigar_op_count(
    dst: &mut Vec<u8>,
    base_count: usize,
    op_count: usize,
    alignment_span: usize,
) -> io::Result<Option<sam::alignment::record_buf::Cigar>> {
    const OVERFLOWING_OP_COUNT: u16 = 2;

    if let Ok(n) = u16::try_from(op_count) {
        write_u16_le(dst, n);
        Ok(None)
    } else {
        write_u16_le(dst, OVERFLOWING_OP_COUNT);

        let k = base_count;
        let m = alignment_span;

        Ok(Some(
            [Op::new(Kind::SoftClip, k), Op::new(Kind::Skip, m)]
                .into_iter()
                .collect(),
        ))
    }
}

pub(super) fn write_cigar(dst: &mut Vec<u8>, cigar: CigarRef<'_>) -> Result<(), EncodeError> {
    match cigar {
        CigarRef::FourBytePacked(src) => write_four_byte_packed_cigar(dst, src),
        CigarRef::Cigar(c) => write_generic_cigar(dst, &c),
    }
}

fn write_four_byte_packed_cigar(dst: &mut Vec<u8>, src: &[u8]) -> Result<(), EncodeError> {
    const MAX_KIND_VALUE: u8 = 8;

    if let (chunks, []) = src.as_chunks() {
        if !chunks
            .iter()
            .all(|[b, _, _, _]| (b & 0x0f) <= MAX_KIND_VALUE)
        {
            return Err(EncodeError::Io(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid kind",
            )));
        }
    } else {
        unreachable!();
    }

    dst.extend(src);

    Ok(())
}

pub(super) fn write_generic_cigar<C>(dst: &mut Vec<u8>, cigar: &C) -> Result<(), EncodeError>
where
    C: Cigar,
{
    for result in cigar.iter() {
        let op = result.map_err(EncodeError::Io)?;
        let n = encode_op(op).map_err(EncodeError::InvalidOp)?;
        write_u32_le(dst, n);
    }

    Ok(())
}

#[cfg(test)]
mod tests {

    use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

    use super::*;

    #[test]
    fn test_overflowing_write_cigar_op_count() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let overflowing_cigar = overflowing_write_cigar_op_count(&mut buf, 4, 1, 4)?;
        assert_eq!(buf, [0x01, 0x00]);
        assert!(overflowing_cigar.is_none());

        const MAX_CIGAR_OP_COUNT: usize = (1 << 16) - 1;
        buf.clear();
        let overflowing_cigar =
            overflowing_write_cigar_op_count(&mut buf, 4, MAX_CIGAR_OP_COUNT + 1, 1 << 16)?;
        assert_eq!(buf, [0x02, 0x00]);
        assert_eq!(
            overflowing_cigar,
            Some(
                [Op::new(Kind::SoftClip, 4), Op::new(Kind::Skip, 1 << 16)]
                    .into_iter()
                    .collect()
            )
        );

        Ok(())
    }

    #[test]
    fn test_write_cigar() -> Result<(), EncodeError> {
        fn t(buf: &mut Vec<u8>, cigar: &CigarBuf, expected: &[u8]) -> Result<(), EncodeError> {
            buf.clear();
            let c = CigarRef::Cigar(Box::new(cigar));
            write_cigar(buf, c)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &CigarBuf::default(), &[])?;
        t(
            &mut buf,
            &[Op::new(Kind::Match, 4)].into_iter().collect(),
            &[0x40, 0x00, 0x00, 0x00],
        )?;
        t(
            &mut buf,
            &[Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)]
                .into_iter()
                .collect(),
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
        )?;

        Ok(())
    }

    #[test]
    fn test_write_four_byte_packed_cigar() -> Result<(), EncodeError> {
        let mut dst = Vec::new();

        dst.clear();
        let src = [0x40, 0x00, 0x00, 0x00];
        write_four_byte_packed_cigar(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00];
        write_four_byte_packed_cigar(&mut dst, &src)?;
        assert_eq!(dst, src);

        dst.clear();
        let src = [0x0f, 0x00, 0x00, 0x00];
        assert!(matches!(
            write_four_byte_packed_cigar(&mut dst, &src),
            Err(EncodeError::Io(e)) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
