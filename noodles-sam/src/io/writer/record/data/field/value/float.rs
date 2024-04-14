use std::io::{self, Write};

use crate::io::writer::num;

pub(super) fn write_float<W>(writer: &mut W, n: f32) -> io::Result<()>
where
    W: Write,
{
    if is_valid(n) {
        num::write_f32(writer, n)
    } else {
        Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid float"))
    }
}

// § 1.5 "The alignment section: optional fields" (2023-05-24): "`[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?`".
fn is_valid(n: f32) -> bool {
    n.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_float() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_float(&mut buf, 0.0)?;
        assert_eq!(buf, b"0");

        buf.clear();
        assert!(matches!(
            write_float(&mut buf, f32::NAN),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_is_valid() {
        assert!(is_valid(0.0));

        assert!(!is_valid(f32::INFINITY));
        assert!(!is_valid(f32::NEG_INFINITY));
        assert!(!is_valid(f32::NAN));
    }
}
