mod op;

use std::{error, fmt, io, mem};

use bytes::Buf;
use noodles_sam::{self as sam, alignment::Record, record::Cigar};

use self::op::decode_op;

/// An error when a raw BAM record CIGAR fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// Unexpected EOF.
    UnexpectedEof,
    /// An op is invalid.
    InvalidOp(op::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::UnexpectedEof => None,
            Self::InvalidOp(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
            Self::InvalidOp(_) => write!(f, "invalid op"),
        }
    }
}

pub(crate) fn get_op_count<B>(src: &mut B) -> Result<usize, ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u16>() {
        return Err(ParseError::UnexpectedEof);
    }

    Ok(usize::from(src.get_u16_le()))
}

pub fn get_cigar<B>(src: &mut B, cigar: &mut Cigar, n_cigar_op: usize) -> Result<(), ParseError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u32>() * n_cigar_op {
        return Err(ParseError::UnexpectedEof);
    }

    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = decode_op(src.get_u32_le()).map_err(ParseError::InvalidOp)?;
        cigar.as_mut().push(op);
    }

    Ok(())
}

// § 4.2.2 "`N_CIGAR_OP` field" (2022-08-22)
pub(super) fn resolve(header: &sam::Header, record: &mut Record) -> io::Result<()> {
    use sam::record::{
        cigar::{op::Kind, Op},
        data::field::{tag, value::Array},
    };

    if let [op_0, op_1] = record.cigar().as_ref() {
        if let Some((_, reference_sequence)) = record.reference_sequence(header).transpose()? {
            let k = record.sequence().len();
            let m = reference_sequence.length().get();

            if *op_0 == Op::new(Kind::SoftClip, k) && *op_1 == Op::new(Kind::Skip, m) {
                if let Some((_, value)) = record.data_mut().remove(&tag::CIGAR) {
                    let data = value
                        .as_array()
                        .and_then(|array| match array {
                            Array::UInt32(values) => Some(values),
                            _ => None,
                        })
                        .ok_or_else(|| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                "invalid CG data field value type",
                            )
                        })?;

                    let cigar = record.cigar_mut();
                    cigar.clear();

                    for &n in data {
                        let op = decode_op(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

                        cigar.as_mut().push(op);
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_op_count() {
        let mut src = &0u16.to_le_bytes()[..];
        assert_eq!(get_op_count(&mut src), Ok(0));

        let mut src = &8u16.to_le_bytes()[..];
        assert_eq!(get_op_count(&mut src), Ok(8));

        let mut src = &[][..];
        assert_eq!(get_op_count(&mut src), Err(ParseError::UnexpectedEof));
    }

    #[test]
    fn test_get_cigar() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            mut src: &[u8],
            actual: &mut Cigar,
            n_cigar_op: usize,
            expected: &Cigar,
        ) -> Result<(), ParseError> {
            get_cigar(&mut src, actual, n_cigar_op)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        let mut buf = Cigar::default();

        t(&[], &mut buf, 0, &Cigar::default())?;
        t(&[0x40, 0x00, 0x00, 0x00], &mut buf, 1, &"4M".parse()?)?;
        t(
            &[0x40, 0x00, 0x00, 0x00, 0x25, 0x00, 0x00, 0x00],
            &mut buf,
            2,
            &"4M2H".parse()?,
        )?;

        Ok(())
    }

    #[test]
    fn test_resolve() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use sam::{
            header::record::value::{map::ReferenceSequence, Map},
            record::{
                cigar::op::{self, Op},
                data::field::{tag, value::Array, Value},
                Cigar,
            },
        };

        let header = sam::Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .build();

        let mut record = Record::builder()
            .set_reference_sequence_id(0)
            .set_cigar("4S8N".parse()?)
            .set_sequence("ACGT".parse()?)
            .set_data(
                [(tag::CIGAR, Value::Array(Array::UInt32(vec![0x40])))]
                    .into_iter()
                    .collect(),
            )
            .build();

        resolve(&header, &mut record)?;

        let expected = Cigar::try_from(vec![Op::new(op::Kind::Match, 4)])?;

        assert_eq!(record.cigar(), &expected);
        assert!(record.data().get(&tag::CIGAR).is_none());

        Ok(())
    }
}
