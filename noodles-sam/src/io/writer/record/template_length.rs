use std::io::{self, Write};

use crate::io::writer::num;

pub(super) fn write_template_length<W>(writer: &mut W, template_length: i32) -> io::Result<()>
where
    W: Write,
{
    if is_valid(template_length) {
        num::write_i32(writer, template_length)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid template length",
        ))
    }
}

// § 1.4 "The alignment section: mandatory fields" (2025-08-12): "[-2^31 + 1, 2^31 - 1]".
fn is_valid(n: i32) -> bool {
    const MIN: i32 = i32::MIN + 1;
    n >= MIN
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_template_length() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_template_length(&mut buf, 0)?;
        assert_eq!(buf, b"0");

        buf.clear();
        assert!(matches!(
            write_template_length(&mut buf, i32::MIN),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(i32::MIN + 1));
        assert!(is_valid(0));
        assert!(is_valid(i32::MAX));

        assert!(!is_valid(i32::MIN));
    }
}
