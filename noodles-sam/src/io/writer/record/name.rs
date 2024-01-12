use std::io::{self, Write};

use crate::alignment::record::Name;

const MAX_LENGTH: usize = 254;
const MISSING: &[u8] = b"*";

pub(super) fn write_name<W, N>(writer: &mut W, name: Option<N>) -> io::Result<()>
where
    W: Write,
    N: Name,
{
    if let Some(name) = name {
        let name = name.as_bytes();

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
    use super::*;
    use crate::alignment::record_buf::Name as NameBuf;

    #[test]
    fn test_write_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, name: Option<&NameBuf>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_name(buf, name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, b"*")?;
        t(&mut buf, Some(&NameBuf::from(b"r0")), b"r0")?;

        Ok(())
    }

    #[test]
    fn test_write_name_with_invalid_name() {
        fn t(raw_name: &[u8]) {
            let mut buf = Vec::new();

            assert!(matches!(
                write_name(&mut buf, Some(&NameBuf::from(raw_name))),
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
