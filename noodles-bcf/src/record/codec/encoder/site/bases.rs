use std::{
    io::{self, Write},
    str,
};

use noodles_vcf::variant::record::{AlternateBases, ReferenceBases};

use crate::record::codec::{Value, encoder::value::write_value};

pub(super) fn write_bases<W, R, A>(
    writer: &mut W,
    reference_bases: R,
    alternate_bases: A,
) -> io::Result<()>
where
    W: Write,
    R: ReferenceBases,
    A: AlternateBases,
{
    write_reference_bases(writer, reference_bases)?;

    if !alternate_bases.is_empty() {
        write_alternate_bases(writer, alternate_bases)?;
    }

    Ok(())
}

fn write_reference_bases<W, R>(writer: &mut W, reference_bases: R) -> io::Result<()>
where
    W: Write,
    R: ReferenceBases,
{
    let raw_bases: Vec<_> = reference_bases
        .iter()
        .map(|result| {
            let b = result?;

            resolve_reference_base(b).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference base: {}", char::from(b)),
                )
            })
        })
        .collect::<io::Result<_>>()?;

    let r#ref =
        str::from_utf8(&raw_bases).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    let ref_value = Some(Value::String(Some(r#ref)));
    write_value(writer, ref_value)?;

    Ok(())
}

// § 1.6.1.4 "Fixed fields: REF" (2026-02-25): "Each base must be one of A,C,G,T,N (case
// insensitive). [...] If the reference sequence contains IUPAC ambiguity codes not allowed by this
// specification (such as R = A/G), the ambiguous reference base must be reduced to a concrete base
// by using the one that is first alphabetically (thus R as a reference base is converted to A in
// VCF.)"
fn resolve_reference_base(b: u8) -> Option<u8> {
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

fn write_alternate_bases<W, A>(writer: &mut W, alternate_bases: A) -> io::Result<()>
where
    W: Write,
    A: AlternateBases,
{
    for result in alternate_bases.iter() {
        let alt = result?;
        let alt_value = Some(Value::String(Some(alt)));
        write_value(writer, alt_value)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_bases() -> Result<(), Box<dyn std::error::Error>> {
        use noodles_vcf::variant::record_buf::AlternateBases as AlternateBasesBuf;

        fn t(
            buf: &mut Vec<u8>,
            reference_bases: &str,
            alternate_bases: &AlternateBasesBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_bases(buf, reference_bases, alternate_bases)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, "A", &AlternateBasesBuf::default(), &[0x17, b'A'])?;
        t(&mut buf, "R", &AlternateBasesBuf::default(), &[0x17, b'A'])?;

        t(
            &mut buf,
            "A",
            &AlternateBasesBuf::from(vec![String::from("G")]),
            &[0x17, b'A', 0x17, b'G'],
        )?;

        t(
            &mut buf,
            "A",
            &AlternateBasesBuf::from(vec![String::from("G"), String::from("T")]),
            &[0x17, b'A', 0x17, b'G', 0x17, b'T'],
        )?;

        Ok(())
    }

    #[test]
    fn test_resolve_reference_base() {
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
            assert_eq!(resolve_reference_base(base), Some(resolved_base));
            assert_eq!(
                resolve_reference_base(base.to_ascii_lowercase()),
                Some(resolved_base.to_ascii_lowercase())
            );
        }

        assert!(resolve_reference_base(b'Z').is_none());
    }
}
