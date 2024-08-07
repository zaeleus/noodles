use std::io::{self, Write};

use bstr::BStr;

const MAX_LENGTH: usize = 254;
const MISSING: &[u8] = b"*";

pub(super) fn write_name<W>(writer: &mut W, name: Option<&BStr>) -> io::Result<()>
where
    W: Write,
{
    if let Some(name) = name {
        if !is_valid(name) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        writer.write_all(name)?;
    } else {
        writer.write_all(MISSING)?;
    }

    Ok(())
}

fn is_valid(buf: &[u8]) -> bool {
    (1..=MAX_LENGTH).contains(&buf.len())
        && buf != MISSING
        && buf.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;

    use super::*;

    #[test]
    fn test_write_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, name: Option<&[u8]>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_name(buf, name.map(|buf| buf.as_bstr()))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b"*")?;
        t(&mut buf, Some(b"r0"), b"r0")?;

        Ok(())
    }

    #[test]
    fn test_write_name_with_invalid_name() {
        fn t(raw_name: &[u8]) {
            let mut buf = Vec::new();

            assert!(matches!(
                write_name(&mut buf, Some(raw_name.as_bstr())),
                Err(e) if e.kind() == io::ErrorKind::InvalidInput
            ));
        }

        t(b"");
        t(b"*");
        t(b"r 0");
        t(b"@r0");

        let s = vec![b'n'; MAX_LENGTH + 1];
        t(&s);
    }
}
