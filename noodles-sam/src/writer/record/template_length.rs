use std::io::{self, Write};

use crate::writer::num;

pub(super) fn write_template_length<W>(writer: &mut W, template_length: i32) -> io::Result<()>
where
    W: Write,
{
    num::write_i32(writer, template_length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_template_length() -> io::Result<()> {
        let mut buf = Vec::new();
        write_template_length(&mut buf, 0)?;
        assert_eq!(buf, b"0");
        Ok(())
    }
}
