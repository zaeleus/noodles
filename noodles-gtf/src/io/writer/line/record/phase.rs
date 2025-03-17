use std::io::{self, Write};

use noodles_gff::feature::record::Phase;

use super::write_missing;

pub(super) fn write_phase<W>(writer: &mut W, phase: Option<Phase>) -> io::Result<()>
where
    W: Write,
{
    const ZERO: &[u8] = b"0";
    const ONE: &[u8] = b"1";
    const TWO: &[u8] = b"2";

    match phase {
        None => write_missing(writer),
        Some(Phase::Zero) => writer.write_all(ZERO),
        Some(Phase::One) => writer.write_all(ONE),
        Some(Phase::Two) => writer.write_all(TWO),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_phase() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, phase: Option<Phase>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_phase(buf, phase)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b".")?;
        t(&mut buf, Some(Phase::Zero), b"0")?;
        t(&mut buf, Some(Phase::One), b"1")?;
        t(&mut buf, Some(Phase::Two), b"2")?;

        Ok(())
    }
}
