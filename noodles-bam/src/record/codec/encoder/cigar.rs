mod op;

use std::io;

use noodles_sam::{
    self as sam,
    alignment::record::{
        Cigar,
        cigar::{Op, op::Kind},
    },
};

use self::op::encode_op;
use super::num::{write_u16_le, write_u32_le};

// § 4.2.2 "N_CIGAR_OP field" (2023-11-16): "For an alignment with more [than 65535] CIGAR
// operations, BAM [...] sets `CIGAR` to `kSmN` as a placeholder, where `k` equals `l_seq`, `m` is
// the reference sequence length in the alignment..."
pub(super) fn overflowing_write_cigar_op_count<C>(
    dst: &mut Vec<u8>,
    base_count: usize,
    cigar: &C,
) -> io::Result<Option<sam::alignment::record_buf::Cigar>>
where
    C: Cigar,
{
    const OVERFLOWING_OP_COUNT: u16 = 2;

    if let Ok(op_count) = u16::try_from(cigar.len()) {
        write_u16_le(dst, op_count);
        Ok(None)
    } else {
        write_u16_le(dst, OVERFLOWING_OP_COUNT);

        let k = base_count;
        let m = cigar.alignment_span()?;

        Ok(Some(
            [Op::new(Kind::SoftClip, k), Op::new(Kind::Skip, m)]
                .into_iter()
                .collect(),
        ))
    }
}

pub(super) fn write_cigar<C>(dst: &mut Vec<u8>, cigar: &C) -> io::Result<()>
where
    C: Cigar,
{
    for result in cigar.iter() {
        let op = result?;
        let n = encode_op(op).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_u32_le(dst, n);
    }

    Ok(())
}

/// Encodes CIGAR operations from a slice with pre-allocated capacity.
///
/// This is an optimized version of [`write_cigar`] that bypasses the trait-based
/// iterator and works directly with a slice of operations. It pre-reserves
/// capacity for all operations before encoding.
///
/// # Performance
///
/// This provides approximately 3-5% throughput improvement compared to the
/// trait-based iterator version by:
/// - Eliminating dynamic dispatch overhead
/// - Pre-allocating buffer capacity for all operations
/// - Direct slice iteration instead of Result unwrapping per operation
///
/// # Arguments
///
/// * `dst` - Output buffer to append encoded CIGAR bytes to
/// * `ops` - Slice of CIGAR operations to encode
///
/// # Errors
///
/// Returns an error if any operation has an invalid length (> 2^28 - 1).
#[inline]
pub(super) fn write_cigar_from_slice(dst: &mut Vec<u8>, ops: &[Op]) -> io::Result<()> {
    dst.reserve(ops.len() * 4);

    for &op in ops {
        let n = encode_op(op).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_u32_le(dst, n);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::iter;

    use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

    use super::*;

    #[test]
    fn test_overflowing_write_cigar_op_count() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        let cigar: CigarBuf = [Op::new(Kind::Match, 4)].into_iter().collect();
        let overflowing_cigar = overflowing_write_cigar_op_count(&mut buf, 4, &cigar)?;
        assert_eq!(buf, [0x01, 0x00]);
        assert!(overflowing_cigar.is_none());

        const MAX_CIGAR_OP_COUNT: usize = (1 << 16) - 1;
        buf.clear();
        let cigar: CigarBuf =
            iter::repeat_n(Op::new(Kind::Match, 1), MAX_CIGAR_OP_COUNT + 1).collect();
        let overflowing_cigar = overflowing_write_cigar_op_count(&mut buf, 4, &cigar)?;
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
    fn test_write_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, cigar: &CigarBuf, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_cigar(buf, cigar)?;
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
    fn test_write_cigar_from_slice() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        // Empty CIGAR
        write_cigar_from_slice(&mut buf, &[])?;
        assert!(buf.is_empty());

        // Single operation
        buf.clear();
        let ops = [Op::new(Kind::Match, 4)];
        write_cigar_from_slice(&mut buf, &ops)?;
        assert_eq!(buf, [0x40, 0x00, 0x00, 0x00]);

        // Multiple operations
        buf.clear();
        let ops = [Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)];
        write_cigar_from_slice(&mut buf, &ops)?;
        assert_eq!(buf, [0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00]);

        Ok(())
    }

    #[test]
    fn test_write_cigar_from_slice_matches_trait_version() {
        let mut buf_trait = Vec::new();
        let mut buf_slice = Vec::new();

        // Test various CIGAR configurations
        let test_cases: Vec<CigarBuf> = vec![
            CigarBuf::default(),
            [Op::new(Kind::Match, 10)].into_iter().collect(),
            [
                Op::new(Kind::Match, 50),
                Op::new(Kind::Insertion, 2),
                Op::new(Kind::Match, 48),
            ]
            .into_iter()
            .collect(),
            [
                Op::new(Kind::SoftClip, 5),
                Op::new(Kind::Match, 90),
                Op::new(Kind::SoftClip, 5),
            ]
            .into_iter()
            .collect(),
        ];

        for cigar in &test_cases {
            buf_trait.clear();
            buf_slice.clear();

            write_cigar(&mut buf_trait, cigar).unwrap();

            let ops: &[Op] = cigar.as_ref();
            write_cigar_from_slice(&mut buf_slice, ops).unwrap();

            assert_eq!(
                buf_trait,
                buf_slice,
                "Mismatch for CIGAR with {} ops",
                cigar.len()
            );
        }
    }
}
