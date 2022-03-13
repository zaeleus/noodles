use std::io::{self, Write};

use noodles_sam::record::data::field::Tag;

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
    fn test_write_tag() -> io::Result<()> {
        let mut buf = Vec::new();
        write_tag(&mut buf, Tag::AlignmentHitCount)?;
        assert_eq!(buf, [b'N', b'H']);
        Ok(())
    }
}
