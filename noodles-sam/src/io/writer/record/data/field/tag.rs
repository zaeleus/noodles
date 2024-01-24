use std::io::{self, Write};

use crate::alignment::record::data::field::Tag;

pub fn write_tag<W>(writer: &mut W, tag: Tag) -> io::Result<()>
where
    W: Write,
{
    if is_valid(tag) {
        writer.write_all(tag.as_ref())
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid tag"))
    }
}

// § 1.5 "The alignment section: optional fields" (2023-05-24): "`[A-Za-z][A-Za-z0-9]`".
fn is_valid(tag: Tag) -> bool {
    let b: [u8; 2] = tag.into();
    b[0].is_ascii_alphabetic() && b[1].is_ascii_alphanumeric()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_tag() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_tag(&mut buf, Tag::ALIGNMENT_HIT_COUNT)?;
        assert_eq!(buf, b"NH");

        buf.clear();
        assert!(matches!(
            write_tag(&mut buf, Tag::new(b'n', b'\t')),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(Tag::ALIGNMENT_HIT_COUNT));
        assert!(is_valid(Tag::new(b'n', b'd')));

        assert!(!is_valid(Tag::new(b'n', b'\t')));
        assert!(!is_valid(Tag::new(b' ', b'd')));
        assert!(!is_valid(Tag::new(0x00, 0xe6)));
    }
}
