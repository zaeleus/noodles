use std::{
    error, fmt,
    io::{self, Write},
};

use crate::variant::record::ReferenceBases;

/// An error returns when reference bases fail to write.
#[derive(Debug)]
pub enum WriteError {
    /// I/O error.
    Io(io::Error),
    /// A reference base is invalid.
    InvalidReferenceBase(u8),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceBase(_) => None,
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceBase(b) => {
                write!(f, "invalid reference base: {}", char::from(*b))
            }
        }
    }
}

pub(super) fn write_reference_bases<W, B>(
    writer: &mut W,
    reference_bases: B,
) -> Result<(), WriteError>
where
    W: Write,
    B: ReferenceBases,
{
    for result in reference_bases.iter() {
        let base = result
            .map_err(WriteError::Io)
            .and_then(|b| resolve_base(b).ok_or(WriteError::InvalidReferenceBase(b)))?;

        writer.write_all(&[base]).map_err(WriteError::Io)?;
    }

    Ok(())
}

// § 1.6.1.4 "Fixed fields: REF" (2024-04-20): "Each base must be one of A,C,G,T,N (case
// insensitive). [...] If the reference sequence contains IUPAC ambiguity codes not allowed by this
// specification (such as R = A/G), the ambiguous reference base must be reduced to a concrete base
// by using the one that is first alphabetically (thus R as a reference base is converted to A in
// VCF.)"
fn resolve_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'W' | b'M' | b'R' | b'D' | b'H' | b'V' => Some(b'A'),
        b'C' | b'S' | b'Y' | b'B' => Some(b'C'),
        b'G' | b'K' => Some(b'G'),
        b'T' => Some(b'T'),
        b'N' => Some(b'N'),

        b'a' | b'w' | b'm' | b'r' | b'd' | b'h' | b'v' => Some(b'a'),
        b'c' | b's' | b'y' | b'b' => Some(b'c'),
        b'g' | b'k' => Some(b'g'),
        b't' => Some(b't'),
        b'n' => Some(b'n'),

        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_bases() -> Result<(), WriteError> {
        let mut buf = Vec::new();

        buf.clear();
        write_reference_bases(&mut buf, "A")?;
        assert_eq!(buf, b"A");

        buf.clear();
        write_reference_bases(&mut buf, "R")?;
        assert_eq!(buf, b"A");

        buf.clear();
        write_reference_bases(&mut buf, "AC")?;
        assert_eq!(buf, b"AC");

        buf.clear();
        assert!(matches!(
            write_reference_bases(&mut buf, "Z"),
            Err(WriteError::InvalidReferenceBase(b'Z'))
        ));

        Ok(())
    }

    #[test]
    fn test_resolve_base() {
        let bases = [
            (b'A', b'A'),
            (b'C', b'C'),
            (b'G', b'G'),
            (b'T', b'T'),
            (b'W', b'A'),
            (b'S', b'C'),
            (b'M', b'A'),
            (b'K', b'G'),
            (b'R', b'A'),
            (b'Y', b'C'),
            (b'B', b'C'),
            (b'D', b'A'),
            (b'H', b'A'),
            (b'V', b'A'),
            (b'N', b'N'),
        ];

        for &(base, resolved_base) in &bases {
            assert_eq!(resolve_base(base), Some(resolved_base));
            assert_eq!(
                resolve_base(base.to_ascii_lowercase()),
                Some(resolved_base.to_ascii_lowercase())
            );
        }

        assert!(resolve_base(b'Z').is_none());
    }
}
