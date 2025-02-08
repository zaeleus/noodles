use std::io::{self, Write};

use byteorder::WriteBytesExt;

use crate::container::block::ContentType;

pub(super) fn write_content_type<W>(writer: &mut W, content_type: ContentType) -> io::Result<()>
where
    W: Write,
{
    let n = match content_type {
        ContentType::FileHeader => 0,
        ContentType::CompressionHeader => 1,
        ContentType::SliceHeader => 2,
        ContentType::Reserved => 3,
        ContentType::ExternalData => 4,
        ContentType::CoreData => 5,
    };

    writer.write_u8(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_content_type() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, content_type: ContentType, expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_content_type(buf, content_type)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, ContentType::FileHeader, &[0x00])?;
        t(&mut buf, ContentType::CompressionHeader, &[0x01])?;
        t(&mut buf, ContentType::SliceHeader, &[0x02])?;
        t(&mut buf, ContentType::Reserved, &[0x03])?;
        t(&mut buf, ContentType::ExternalData, &[0x04])?;
        t(&mut buf, ContentType::CoreData, &[0x05])?;

        Ok(())
    }
}
