use std::io::{self, Write};

use super::{Header, Record};

/// A SAM writer.
///
/// The SAM format is comprised of two parts: 1) a header and 2) a list of records.
///
/// Each header line is prefixed with an `@` (at sign). The header is optional and may be empty.
///
/// SAM records are line-based and follow directly after the header or the start of the file until
/// EOF.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_sam as sam;
///
/// let mut writer = sam::Writer::new(Vec::new());
///
/// let header = sam::Header::builder().add_comment("noodles-sam").build();
/// writer.write_header(&header)?;
///
/// let record = sam::Record::default();
/// writer.write_record(&record)?;
///
/// let expected = b"@CO\tnoodles-sam
/// *\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*
/// ";
///
/// assert_eq!(&writer.get_ref()[..], &expected[..]);
/// # Ok::<(), io::Error>(())
/// ```
pub struct Writer<W>
where
    W: Write,
{
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a SAM writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let writer = sam::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a SAM header.
    ///
    /// The SAM header is optional, though recommended to include. A call to this method can be
    /// omitted if it is empty.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// let mut writer = sam::Writer::new(Vec::new());
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// writer.write_header(&header)?;
    /// assert_eq!(writer.get_ref(), b"@CO\tnoodles-sam\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write!(self.inner, "{}", header)
    }

    /// Writes a SAM record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_sam as sam;
    /// let mut writer = sam::Writer::new(Vec::new());
    /// let record = sam::Record::default();
    /// writer.write_record(&record)?;
    /// assert_eq!(writer.get_ref(), b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        validate(record)?;
        writeln!(self.inner, "{}", record)
    }
}

/// Validate SAM field relations
pub fn validate(record: &Record) -> io::Result<()> {
    // Before writing validate record requirements for BAM encoding
    // Sequence and quality must be same length
    // If no qualities are present use 0xFF as placeholders
    // https://samtools.github.io/hts-specs/SAMv1.pdf#subsubsection.4.2.3
    let sequence_len = record.sequence().len();
    let quality_scores_len = record.quality_scores().len();

    if quality_scores_len > 0 && quality_scores_len != sequence_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "quality scores length mismatch: expected {}, got {}",
                sequence_len, quality_scores_len
            ),
        ));
    }

    // The sum of the query consuming CIGAR operations must equal the SEQ length
    // https://samtools.github.io/hts-specs/SAMv1.pdf#subsubsection.4.2.2
    let sequence_len =
        u32::try_from(sequence_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    let cigar_query_len = record.cigar().read_len();
    if cigar_query_len != sequence_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "cigar operations length mismatch: expected {}, got {}",
                sequence_len, cigar_query_len
            ),
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::record::{
        data::{self, Field},
        Data,
    };

    use super::*;

    #[test]
    fn test_write_record_with_data() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(vec![]);

        let data = Data::try_from(vec![Field::new(
            data::field::Tag::ReadGroup,
            data::field::Value::String(String::from("rg0")),
        )])?;

        let record = Record::builder().set_data(data).build()?;

        writer.write_record(&record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";

        assert_eq!(&writer.get_ref()[..], &expected[..]);

        Ok(())
    }
}
