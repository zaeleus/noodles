use std::io;

use bytes::BufMut;
use noodles_sam::alignment::record::Name;

const MAX_LENGTH: usize = 254;
const MISSING: &[u8] = b"*";

pub fn put_name<B, N>(dst: &mut B, name: Option<N>) -> io::Result<()>
where
    B: BufMut,
    N: Name,
{
    const NUL: u8 = 0x00;

    if let Some(name) = name {
        let buf = name.as_bytes();

        if !is_valid(buf) {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }

        dst.put(buf);
    } else {
        dst.put(MISSING);
    }

    dst.put_u8(NUL);

    Ok(())
}

fn is_valid(buf: &[u8]) -> bool {
    (1..=MAX_LENGTH).contains(&buf.len())
        && buf != MISSING
        && buf.iter().all(|&b| b.is_ascii_graphic() && b != b'@')
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::record_buf::Name as NameBuf;

    use super::*;

    #[test]
    fn test_put_name() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, name: Option<&NameBuf>, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            put_name(buf, name)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, None, &[b'*', 0x00])?;
        t(&mut buf, Some(&NameBuf::from(b"r0")), &[b'r', b'0', 0x00])?;

        Ok(())
    }

    #[test]
    fn test_put_name_with_invalid_name() {
        fn t(raw_name: &[u8]) {
            let mut buf = Vec::new();
            let name = NameBuf::from(raw_name);
            assert!(matches!(
                put_name(&mut buf, Some(&name)),
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
