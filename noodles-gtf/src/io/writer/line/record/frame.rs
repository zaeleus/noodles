use std::io::{self, Write};

use super::write_missing;
use crate::record_buf::Frame;

pub(super) fn write_frame<W>(writer: &mut W, frame: Option<Frame>) -> io::Result<()>
where
    W: Write,
{
    if let Some(frame) = frame {
        write!(writer, "{frame}")
    } else {
        write_missing(writer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_frame() -> Result<(), Box<dyn std::error::Error>> {
        fn t(buf: &mut Vec<u8>, frame: Option<Frame>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_frame(buf, frame)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, Some(Frame::try_from(0)?), b"0")?;
        t(&mut buf, Some(Frame::try_from(1)?), b"1")?;
        t(&mut buf, Some(Frame::try_from(2)?), b"2")?;
        t(&mut buf, None, b".")?;

        Ok(())
    }
}
