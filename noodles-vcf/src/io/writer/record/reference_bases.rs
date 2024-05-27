use std::io::{self, Write};

use crate::variant::record::ReferenceBases;

pub(super) fn write_reference_bases<W, B>(writer: &mut W, reference_bases: B) -> io::Result<()>
where
    W: Write,
    B: ReferenceBases,
{
    for result in reference_bases.iter() {
        let base = result.and_then(resolve_base)?;
        writer.write_all(&[base])?;
    }

    Ok(())
}

// § 1.6.1.4 "Fixed fields: REF" (2024-04-20): "Each base must be one of A,C,G,T,N (case
// insensitive). [...] If the reference sequence contains IUPAC ambiguity codes not allowed by this
// specification (such as R = A/G), the ambiguous reference base must be reduced to a concrete base
// by using the one that is first alphabetically (thus R as a reference base is converted to A in
// VCF.)"
fn resolve_base(b: u8) -> io::Result<u8> {
    match b {
        b'A' | b'W' | b'M' | b'R' | b'D' | b'H' | b'V' => Ok(b'A'),
        b'C' | b'S' | b'Y' | b'B' => Ok(b'C'),
        b'G' | b'K' => Ok(b'G'),
        b'T' => Ok(b'T'),
        b'N' => Ok(b'N'),

        b'a' | b'w' | b'm' | b'r' | b'd' | b'h' | b'v' => Ok(b'a'),
        b'c' | b's' | b'y' | b'b' => Ok(b'c'),
        b'g' | b'k' => Ok(b'g'),
        b't' => Ok(b't'),
        b'n' => Ok(b'n'),

        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid reference base",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_bases() -> io::Result<()> {
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
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_resolve_base() -> io::Result<()> {
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
            assert_eq!(resolve_base(base)?, resolved_base);
            assert_eq!(
                resolve_base(base.to_ascii_lowercase())?,
                resolved_base.to_ascii_lowercase()
            );
        }

        Ok(())
    }
}
