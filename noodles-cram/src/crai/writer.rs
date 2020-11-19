use std::io::{self, Write};

use flate2::write::GzEncoder;
use noodles_bam as bam;

use super::Record;

/// A CRAM index writer.
pub struct Writer<W>
where
    W: Write,
{
    inner: GzEncoder<W>,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a CRAM index writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let writer = crai::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self {
            inner: GzEncoder::new(inner, Default::default()),
        }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::crai;
    /// let writer = crai::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        self.inner.get_ref()
    }

    /// Attempts to finish the output stream and returns the underlying writer.
    ///
    /// This is typically only manually called if the underlying stream is needed before the writer
    /// is dropped.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::crai;
    /// let writer = crai::Writer::new(Vec::new());
    /// let empty_gz = [31, 139, 8, 0, 0, 0, 0, 0, 0, 255, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    /// assert_eq!(writer.finish()?, empty_gz);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn finish(self) -> io::Result<W> {
        self.inner.finish()
    }

    /// Writes a CRAM index record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
    /// use noodles_bam as bam;
    /// use noodles_cram::crai;
    ///
    /// let mut writer = crai::Writer::new(Vec::new());
    ///
    /// let record = crai::Record::new(
    ///     bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
    ///     10946,
    ///     6765,
    ///     17711,
    ///     233,
    ///     317811,
    /// );
    ///
    /// writer.write_record(&record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        let reference_sequence_id = record
            .reference_sequence_id()
            .map(i32::from)
            .unwrap_or(bam::record::reference_sequence_id::UNMAPPED);

        writeln!(
            self.inner,
            "{}\t{}\t{}\t{}\t{}\t{}",
            reference_sequence_id,
            record.alignment_start(),
            record.alignment_span(),
            record.offset(),
            record.landmark(),
            record.slice_length()
        )
    }
}

#[cfg(test)]
mod tests {
    use std::{convert::TryFrom, io::Read};

    use crate::crai;

    use super::*;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let record = crai::Record::new(
            bam::record::ReferenceSequenceId::try_from(0).map(Some)?,
            10946,
            6765,
            17711,
            233,
            317811,
        );

        writer.write_record(&record)?;

        let data = writer.finish()?;
        let mut decoder = flate2::read::GzDecoder::new(&data[..]);

        let mut uncompressed_data = String::new();
        decoder.read_to_string(&mut uncompressed_data)?;

        let expected = "0\t10946\t6765\t17711\t233\t317811\n";

        assert_eq!(uncompressed_data, expected);

        Ok(())
    }
}
