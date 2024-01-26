use std::io::{self, Write};

use crate::header::FileFormat;

pub(crate) fn write_file_format<W>(writer: &mut W, file_format: FileFormat) -> io::Result<()>
where
    W: Write,
{
    const PREFIX: &str = "VCFv";
    const DELIMITER: char = '.';

    write!(
        writer,
        "{}{}{}{}",
        PREFIX,
        file_format.major(),
        DELIMITER,
        file_format.minor()
    )?;

    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_write_file_format() -> io::Result<()> {
        let mut buf = Vec::new();
        let file_format = FileFormat::new(4, 3);
        write_file_format(&mut buf, file_format)?;
        assert_eq!(buf, b"VCFv4.3");
        Ok(())
    }
}
