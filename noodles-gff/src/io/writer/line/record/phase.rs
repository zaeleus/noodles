use std::io::{self, Write};

use super::write_missing;
use crate::record_buf::Phase;

const ZERO: &[u8] = b"0";
const ONE: &[u8] = b"1";
const TWO: &[u8] = b"2";

pub(super) fn write_phase<W>(writer: &mut W, phase: Option<Phase>) -> io::Result<()>
where
    W: Write,
{
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
    fn test_write_phase() -> io::Result<()> {
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
