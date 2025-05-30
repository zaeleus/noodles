mod op;

use std::{error, fmt, mem};

use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        record::cigar::{Op, op::Kind},
        record_buf::Cigar,
    },
};

pub(crate) use self::op::decode_op;

/// An error when a raw BAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// An op is invalid.
    InvalidOp(op::DecodeError),
    /// The `CG` data field type is invalid.
    InvalidDataType,
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidOp(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
            Self::InvalidDataType => write!(f, "invalid CG data field type"),
        }
    }
}

pub(super) fn read_op_count(src: &mut &[u8]) -> Result<usize, DecodeError> {
    read_u16_le(src).map(usize::from)
}

pub(super) fn read_cigar(
    src: &mut &[u8],
    cigar: &mut Cigar,
    op_count: usize,
) -> Result<(), DecodeError> {
    let len = mem::size_of::<u32>() * op_count;
    let (mut buf, rest) = src
        .split_at_checked(len)
        .ok_or(DecodeError::UnexpectedEof)?;

    *src = rest;

    let dst = cigar.as_mut();
    dst.clear();

    for _ in 0..op_count {
        let n = read_u32_le(&mut buf)?;
        let op = decode_op(n).map_err(DecodeError::InvalidOp)?;
        dst.push(op);
    }

    Ok(())
}

// § 4.2.2 "`N_CIGAR_OP` field" (2022-08-22)
pub(super) fn resolve(record: &mut RecordBuf) -> Result<(), DecodeError> {
    use sam::alignment::{
        record::data::field::Tag,
        record_buf::data::field::{Value, value::Array},
    };

    if let [op_0, op_1] = record.cigar().as_ref() {
        let k = record.sequence().len();

        if *op_0 == Op::new(Kind::SoftClip, k) && op_1.kind() == Kind::Skip {
            if let Some((_, value)) = record.data_mut().remove(&Tag::CIGAR) {
                let Value::Array(Array::UInt32(values)) = value else {
                    return Err(DecodeError::InvalidDataType);
                };

                let cigar = record.cigar_mut().as_mut();
                cigar.clear();

                for n in values {
                    let op = decode_op(n).map_err(DecodeError::InvalidOp)?;
                    cigar.push(op);
                }
            }
        }
    }

    Ok(())
}

fn read_u16_le(src: &mut &[u8]) -> Result<u16, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u16::from_le_bytes(*buf))
}

fn read_u32_le(src: &mut &[u8]) -> Result<u32, DecodeError> {
    let (buf, rest) = src.split_first_chunk().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(u32::from_le_bytes(*buf))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_op_count() {
        let mut src = &0u16.to_le_bytes()[..];
        assert_eq!(read_op_count(&mut src), Ok(0));

        let mut src = &8u16.to_le_bytes()[..];
        assert_eq!(read_op_count(&mut src), Ok(8));

        let mut src = &[][..];
        assert_eq!(read_op_count(&mut src), Err(DecodeError::UnexpectedEof));
    }

    #[test]
    fn test_read_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut Cigar,
            n_cigar_op: usize,
            expected: &Cigar,
        ) -> Result<(), DecodeError> {
            read_cigar(&mut src, actual, n_cigar_op)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Cigar::default();

        t(&[], &mut buf, 0, &Cigar::default())?;
        t(
            &[0x40, 0x00, 0x00, 0x00],
            &mut buf,
            1,
            &[Op::new(Kind::Match, 4)].into_iter().collect(),
        )?;
        t(
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
            &mut buf,
            2,
            &[Op::new(Kind::Match, 4), Op::new(Kind::HardClip, 2)]
                .into_iter()
                .collect(),
        )?;

        Ok(())
    }

    #[test]
    fn test_resolve() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::{
            record::{cigar::Op, data::field::Tag},
            record_buf::{
                Sequence,
                data::field::{Value, value::Array},
            },
        };

        let mut record = RecordBuf::builder()
            .set_reference_sequence_id(0)
            .set_cigar(
                [Op::new(Kind::SoftClip, 4), Op::new(Kind::Skip, 8)]
                    .into_iter()
                    .collect(),
            )
            .set_sequence(Sequence::from(b"ACGT"))
            .set_data(
                [(Tag::CIGAR, Value::Array(Array::UInt32(vec![0x40])))]
                    .into_iter()
                    .collect(),
            )
            .build();

        resolve(&mut record)?;

        let expected = [Op::new(Kind::Match, 4)].into_iter().collect();

        assert_eq!(record.cigar(), &expected);
        assert!(record.data().get(&Tag::CIGAR).is_none());

        Ok(())
    }
}
