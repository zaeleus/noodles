use std::io::{self, Write};

pub(super) fn write_key<W, K>(writer: &mut W, key: K) -> io::Result<()>
where
    W: Write,
    K: AsRef<str>,
{
    writer.write_all(key.as_ref().as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_key() -> io::Result<()> {
        let mut buf = Vec::new();
        write_key(&mut buf, "INFO")?;
        assert_eq!(buf, b"INFO");
        Ok(())
    }
}
