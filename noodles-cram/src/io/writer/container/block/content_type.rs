use std::io::{self, Write};

use crate::{container::block::ContentType, io::writer::num::write_u8};

pub(super) fn write_content_type<W>(writer: &mut W, content_type: ContentType) -> io::Result<()>
where
    W: Write,
{
    let n = encode(content_type);
    write_u8(writer, n)
}

fn encode(content_type: ContentType) -> u8 {
    match content_type {
        ContentType::FileHeader => 0,
        ContentType::CompressionHeader => 1,
        ContentType::SliceHeader => 2,
        ContentType::Reserved => 3,
        ContentType::ExternalData => 4,
        ContentType::CoreData => 5,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_content_type() -> io::Result<()> {
        let mut buf = Vec::new();
        write_content_type(&mut buf, ContentType::FileHeader)?;
        assert_eq!(buf, [0x00]);
        Ok(())
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(ContentType::FileHeader), 0);
        assert_eq!(encode(ContentType::CompressionHeader), 1);
        assert_eq!(encode(ContentType::SliceHeader), 2);
        assert_eq!(encode(ContentType::Reserved), 3);
        assert_eq!(encode(ContentType::ExternalData), 4);
        assert_eq!(encode(ContentType::CoreData), 5);
    }
}
