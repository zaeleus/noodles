use std::io;

use crate::container::block::ContentType;

pub(super) fn read_content_type(src: &mut &[u8]) -> io::Result<ContentType> {
    let (n, rest) = src
        .split_first()
        .ok_or_else(|| io::Error::from(io::ErrorKind::UnexpectedEof))?;

    *src = rest;

    decode(*n)
}

pub(crate) fn decode(n: u8) -> io::Result<ContentType> {
    match n {
        0 => Ok(ContentType::FileHeader),
        1 => Ok(ContentType::CompressionHeader),
        2 => Ok(ContentType::SliceHeader),
        3 => Ok(ContentType::Reserved),
        4 => Ok(ContentType::ExternalData),
        5 => Ok(ContentType::CoreData),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid content type",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_content_type() -> io::Result<()> {
        assert_eq!(
            read_content_type(&mut &[0x00][..])?,
            ContentType::FileHeader
        );

        assert!(matches!(
            read_content_type(&mut &[][..]),
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }

    #[test]
    fn test_decode() -> io::Result<()> {
        assert_eq!(decode(0)?, ContentType::FileHeader);
        assert_eq!(decode(1)?, ContentType::CompressionHeader);
        assert_eq!(decode(2)?, ContentType::SliceHeader);
        assert_eq!(decode(3)?, ContentType::Reserved);
        assert_eq!(decode(4)?, ContentType::ExternalData);
        assert_eq!(decode(5)?, ContentType::CoreData);

        assert!(matches!(
            decode(6),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
