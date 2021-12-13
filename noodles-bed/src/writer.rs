use std::{
    fmt,
    io::{self, Write},
};

use super::Record;

/// A BED writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a BED writer.
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a BED record.
    pub fn write_record<const N: u8>(&mut self, record: &Record<N>) -> io::Result<()>
    where
        Record<N>: fmt::Display,
    {
        write_record(&mut self.inner, record)
    }
}

fn write_record<W, const N: u8>(writer: &mut W, record: &Record<N>) -> io::Result<()>
where
    W: Write,
    Record<N>: fmt::Display,
{
    writeln!(writer, "{}", record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let record: Record<3> = "sq0\t8\t13".parse()?;
        write_record(&mut buf, &record)?;
        assert_eq!(buf, b"sq0\t8\t13\n");
        Ok(())
    }
}
