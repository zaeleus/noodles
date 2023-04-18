use std::{io, num::NonZeroUsize};

use bytes::Buf;
use noodles_sam::record::ReadName;

pub fn get_read_name<B>(
    src: &mut B,
    read_name: &mut Option<ReadName>,
    l_read_name: NonZeroUsize,
) -> io::Result<()>
where
    B: Buf,
{
    const NUL: u8 = 0x00;
    const MISSING: [u8; 2] = [b'*', NUL];

    let len = usize::from(l_read_name);

    if src.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    *read_name = if src.take(len).chunk() == MISSING {
        src.advance(MISSING.len());
        None
    } else {
        let mut dst = read_name.take().map(Vec::from).unwrap_or_default();

        // SAFETY: len is guaranteed to be > 0.
        dst.resize(len - 1, 0);
        src.copy_to_slice(&mut dst);

        let terminator = src.get_u8();

        if terminator != NUL {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid read name terminator: expected {NUL:#04x}, got {terminator:#04x}"),
            ));
        }

        ReadName::try_from(dst)
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
        fn t(mut src: &[u8], expected: Option<ReadName>) -> io::Result<()> {
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

        // An invalid NUL-terminator.
        let data = [b'*', b'*'];
        let mut src = &data[..];
        let l_read_name = NonZeroUsize::try_from(data.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        assert!(matches!(
            get_read_name(&mut src, &mut None, l_read_name),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        // An invalid character.
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
