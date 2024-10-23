use std::{error, fmt};

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

pub(super) fn read_mapping_quality(src: &mut &[u8]) -> Result<Option<MappingQuality>, DecodeError> {
    let (n, rest) = src.split_first().ok_or(DecodeError::UnexpectedEof)?;
    *src = rest;
    Ok(MappingQuality::new(*n))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_mapping_quality() {
        fn t(mut buf: &[u8], expected: Option<MappingQuality>) {
            assert_eq!(read_mapping_quality(&mut buf), Ok(expected));
        }

        t(&[0x00], MappingQuality::new(0));
        t(&[0x08], MappingQuality::new(8));
        t(&[0xff], None);

        let mut src = &[][..];
        assert_eq!(
            read_mapping_quality(&mut src),
            Err(DecodeError::UnexpectedEof)
        );
    }
}
