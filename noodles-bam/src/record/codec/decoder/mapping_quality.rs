use std::{error, fmt, mem};

use bytes::Buf;
use noodles_sam::alignment::record::MappingQuality;

/// An error when a raw BAM record mapping quality fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// Unexpected EOF.
    UnexpectedEof,
}

impl error::Error for DecodeError {}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedEof => write!(f, "unexpected EOF"),
        }
    }
}

pub(super) fn get_mapping_quality<B>(src: &mut B) -> Result<Option<MappingQuality>, DecodeError>
where
    B: Buf,
{
    if src.remaining() < mem::size_of::<u8>() {
        return Err(DecodeError::UnexpectedEof);
    }

    Ok(MappingQuality::new(src.get_u8()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_mapping_quality() {
        fn t(mut buf: &[u8], expected: Option<MappingQuality>) {
            assert_eq!(get_mapping_quality(&mut buf), Ok(expected));
        }

        t(&[0x00], MappingQuality::new(0));
        t(&[0x08], MappingQuality::new(8));
        t(&[0xff], None);

        let mut src = &[][..];
        assert_eq!(
            get_mapping_quality(&mut src),
            Err(DecodeError::UnexpectedEof)
        );
    }
}
