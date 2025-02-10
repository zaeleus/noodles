use std::io::{self, Write};

pub(super) fn write_file_id<W>(writer: &mut W, file_id: &[u8]) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(file_id)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_file_id() -> io::Result<()> {
        let mut buf = Vec::new();

        let file_id = vec![0; 20];
        write_file_id(&mut buf, &file_id)?;

        assert_eq!(buf, file_id);

        Ok(())
    }
}
