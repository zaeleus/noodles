use std::io::{self, Write};

use crate::record::data::field::Tag;

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
        use crate::record::data::field::tag;

        let mut buf = Vec::new();
        write_tag(&mut buf, tag::ALIGNMENT_HIT_COUNT)?;
        assert_eq!(buf, b"NH");

        Ok(())
    }
}
