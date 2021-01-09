use std::{
    cmp,
    io::{self, Write},
};

use super::{reader::DEFINITION_PREFIX, Record};

const LINE_BASES: u64 = 80;

/// A FASTA writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a FASTA writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let writer = fasta::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// let writer = fasta::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a FASTA record.
    ///
    /// Sequence lines are hard wrapped at 80 bases.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_fasta as fasta;
    ///
    /// let mut writer = fasta::Writer::new(Vec::new());
    ///
    /// let definition = fasta::record::Definition::new(String::from("sq0"), None);
    /// let record = fasta::Record::new(definition, b"ACGT".to_vec());
    ///
    /// writer.write_record(&record)?;
    ///
    /// assert_eq!(writer.get_ref(), b">sq0\nACGT\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record_definition(
            &mut self.inner,
            record.reference_sequence_name(),
            record.description(),
        )?;

        write_record_sequence(&mut self.inner, record.sequence(), LINE_BASES as usize)?;

        Ok(())
    }
}

fn write_record_definition<W>(
    writer: &mut W,
    reference_sequence_name: &str,
    description: Option<&str>,
) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[DEFINITION_PREFIX])?;
    writer.write_all(reference_sequence_name.as_bytes())?;

    if let Some(d) = description {
        writer.write_all(&[b' '])?;
        writer.write_all(d.as_bytes())?;
    }

    writeln!(writer)
}

fn write_record_sequence<W>(writer: &mut W, sequence: &[u8], line_bases: usize) -> io::Result<()>
where
    W: Write,
{
    let mut start = 0;

    while start < sequence.len() {
        let end = cmp::min(start + line_bases, sequence.len());
        let line = &sequence[start..end];

        writer.write_all(line)?;
        writeln!(writer)?;

        start = end;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record_definition() -> io::Result<()> {
        let mut writer = Vec::new();
        write_record_definition(&mut writer, "sq0", None)?;
        assert_eq!(writer, b">sq0\n");

        writer.clear();
        write_record_definition(&mut writer, "sq0", Some("LN:4"))?;
        assert_eq!(writer, b">sq0 LN:4\n");

        Ok(())
    }

    #[test]
    fn test_write_record_sequence() -> io::Result<()> {
        let mut writer = Vec::new();
        write_record_sequence(&mut writer, b"AC", 4)?;
        assert_eq!(writer, b"AC\n");

        writer.clear();
        write_record_sequence(&mut writer, b"ACGT", 4)?;
        assert_eq!(writer, b"ACGT\n");

        writer.clear();
        write_record_sequence(&mut writer, b"ACGTACGT", 4)?;
        assert_eq!(writer, b"ACGT\nACGT\n");

        writer.clear();
        write_record_sequence(&mut writer, b"ACGTACGTAC", 4)?;
        assert_eq!(writer, b"ACGT\nACGT\nAC\n");

        Ok(())
    }
}
