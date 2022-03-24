use std::{io, num::NonZeroUsize};

use bytes::Buf;
use noodles_sam as sam;

pub(super) fn get_read_name<B>(
    buf: &mut B,
    read_name: &mut Option<sam::record::ReadName>,
    l_read_name: NonZeroUsize,
) -> io::Result<()>
where
    B: Buf,
{
    const MISSING: u8 = b'*';

    let len = usize::from(l_read_name);

    if buf.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    *read_name = if len == 2 && buf.chunk()[0] == MISSING {
        buf.advance(2);
        None
    } else {
        let mut read_name_buf = read_name.take().map(Vec::from).unwrap_or_default();

        // SAFETY: len is guaranteed to be > 0.
        read_name_buf.resize(len - 1, Default::default());
        buf.copy_to_slice(&mut read_name_buf);

        // Discard the NUL terminator.
        buf.advance(1);

        sam::record::ReadName::try_from(read_name_buf)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    };

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_read_name() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], expected: Option<sam::record::ReadName>) -> io::Result<()> {
            let mut actual = None;
            let l_read_name = NonZeroUsize::try_from(src.len())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            get_read_name(&mut src, &mut actual, l_read_name)?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'*', 0x00], None)?;
        t(&[b'r', 0x00], "r".parse().map(Some)?)?;
        t(&[b'r', b'1', 0x00], "r1".parse().map(Some)?)?;

        let data = [0xf0, 0x9f, 0x8d, 0x9c, 0x00]; // "🍜\x00"
        let mut src = &data[..];
        let l_read_name = NonZeroUsize::try_from(data.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        assert!(matches!(
            get_read_name(&mut src, &mut None, l_read_name),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
