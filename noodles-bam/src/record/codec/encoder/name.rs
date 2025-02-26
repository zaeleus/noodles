use std::{io, mem};

use bstr::BStr;

use super::num::write_u8;

const MAX_LENGTH: usize = 254;
pub(super) const MISSING: &[u8] = b"*";

pub(super) fn write_length(dst: &mut Vec<u8>, name: Option<&BStr>) -> io::Result<()> {
    let mut len = name.map(|s| s.len()).unwrap_or(MISSING.len());

    // + NUL terminator
    len += mem::size_of::<u8>();

    let n = u8::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_u8(dst, n);

    Ok(())
}

pub(super) fn write_name(dst: &mut Vec<u8>, name: Option<&BStr>) -> io::Result<()> {
    const NUL: u8 = 0x00;

    if let Some(name) = name {
        if !is_valid(name) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        dst.extend_from_slice(name.as_ref());
    } else {
        dst.extend(MISSING);
    }

    write_u8(dst, NUL);

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
    fn test_write_length() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_length(&mut buf, None)?;
        assert_eq!(buf, [0x02]);

        buf.clear();
        write_length(&mut buf, Some(b"r0".as_bstr()))?;
        assert_eq!(buf, [0x03]);

        buf.clear();
        let name = vec![b'n'; 255];
        assert!(matches!(
            write_length(&mut buf, Some(name.as_bstr())),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, name: Option<&[u8]>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_name(buf, name.map(|buf| buf.as_bstr()))?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[b'*', 0x00])?;
        t(&mut buf, Some(b"r0"), &[b'r', b'0', 0x00])?;

        Ok(())
    }

    #[test]
    fn test_put_name_with_invalid_name() {
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
