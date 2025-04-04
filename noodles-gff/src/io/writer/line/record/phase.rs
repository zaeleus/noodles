use std::io::{self, Write};

use bstr::BStr;

use super::write_missing;
use crate::feature::record::Phase;

const CODING_SEQUENCE: &[u8] = b"CDS";

const ZERO: &[u8] = b"0";
const ONE: &[u8] = b"1";
const TWO: &[u8] = b"2";

pub(super) fn write_phase<W>(writer: &mut W, ty: &BStr, phase: Option<Phase>) -> io::Result<()>
where
    W: Write,
{
    match phase {
        None => {
            if ty == CODING_SEQUENCE {
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "missing phase for CDS",
                ))
            } else {
                write_missing(writer)
            }
        }
        Some(Phase::Zero) => writer.write_all(ZERO),
        Some(Phase::One) => writer.write_all(ONE),
        Some(Phase::Two) => writer.write_all(TWO),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_phase() -> io::Result<()> {
        const MISSING: &[u8] = b".";

        fn t(
            buf: &mut Vec<u8>,
            ty: &BStr,
            phase: Option<Phase>,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_phase(buf, ty, phase)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, BStr::new(MISSING), None, b".")?;
        t(&mut buf, BStr::new(MISSING), Some(Phase::Zero), b"0")?;
        t(&mut buf, BStr::new(MISSING), Some(Phase::One), b"1")?;
        t(&mut buf, BStr::new(MISSING), Some(Phase::Two), b"2")?;

        buf.clear();
        assert!(matches!(
            write_phase(&mut buf, BStr::new(CODING_SEQUENCE), None),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
