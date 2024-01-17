use std::io::{self, Write};

use crate::alignment::record::data::field::Tag;

pub fn write_tag<W>(writer: &mut W, tag: Tag) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(tag.as_ref())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() -> io::Result<()> {
        let mut buf = Vec::new();
        write_tag(&mut buf, Tag::ALIGNMENT_HIT_COUNT)?;
        assert_eq!(buf, b"NH");
        Ok(())
    }
}
