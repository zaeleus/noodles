mod num;
mod record;

use std::io::{self, Write};

use self::record::{write_cigar, write_data, write_position, write_quality_scores, write_sequence};
use super::{alignment::Record, AlignmentWriter, Header};

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
/// use noodles_sam::{self as sam, alignment::Record};
///
/// let mut writer = sam::Writer::new(Vec::new());
///
/// let header = sam::Header::builder().add_comment("noodles-sam").build();
/// writer.write_header(&header)?;
///
/// let record = Record::default();
/// writer.write_record(&header, &record)?;
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

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam as sam;
    /// let mut writer = sam::Writer::new(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
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
    /// use noodles_sam::{self as sam, alignment::Record};
    ///
    /// let mut writer = sam::Writer::new(Vec::new());
    ///
    /// let header = sam::Header::default();
    /// let record = Record::default();
    /// writer.write_record(&header, &record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &Header, record: &Record) -> io::Result<()> {
        self.write_alignment_record(header, record)
    }
}

impl<W> AlignmentWriter for Writer<W>
where
    W: Write,
{
    fn write_alignment_header(&mut self, header: &Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_alignment_record(&mut self, header: &Header, record: &Record) -> io::Result<()> {
        const DELIMITER: &[u8] = b"\t";
        const EQ: &[u8] = b"=";
        const MISSING: &[u8] = b"*";

        let qname = record
            .read_name()
            .map(|name| name.as_ref())
            .unwrap_or(MISSING);

        let reference_sequence = record.reference_sequence(header).transpose()?;

        let rname = reference_sequence
            .map(|rs| rs.name().as_bytes())
            .unwrap_or(MISSING);

        let mapq = record
            .mapping_quality()
            .map(u8::from)
            .unwrap_or(crate::record::mapping_quality::MISSING);

        let rnext = record
            .mate_reference_sequence(header)
            .transpose()?
            .map(|mate_reference_sequence| {
                if let Some(rs) = reference_sequence {
                    if mate_reference_sequence.name() == rs.name() {
                        return EQ;
                    }
                }

                mate_reference_sequence.name().as_ref()
            })
            .unwrap_or(MISSING);

        self.inner.write_all(qname)?;

        self.inner.write_all(DELIMITER)?;
        num::write_u16(&mut self.inner, u16::from(record.flags()))?;

        self.inner.write_all(DELIMITER)?;
        self.inner.write_all(rname)?;

        self.inner.write_all(DELIMITER)?;
        write_position(&mut self.inner, record.alignment_start())?;

        self.inner.write_all(DELIMITER)?;
        num::write_u8(&mut self.inner, mapq)?;

        self.inner.write_all(DELIMITER)?;
        write_cigar(&mut self.inner, record.cigar())?;

        self.inner.write_all(DELIMITER)?;
        self.inner.write_all(rnext)?;

        self.inner.write_all(DELIMITER)?;
        write_position(&mut self.inner, record.mate_alignment_start())?;

        self.inner.write_all(DELIMITER)?;
        num::write_i32(&mut self.inner, record.template_length())?;

        self.inner.write_all(DELIMITER)?;
        write_sequence(&mut self.inner, record.sequence())?;

        self.inner.write_all(DELIMITER)?;
        write_quality_scores(&mut self.inner, record.quality_scores())?;

        write_data(&mut self.inner, record.data())?;

        writeln!(self.inner)?;

        Ok(())
    }

    fn finish(&mut self, _: &Header) -> io::Result<()> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record_with_data() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::{
            data::{
                field::{Tag, Value},
                Field,
            },
            Data,
        };

        let mut writer = Writer::new(vec![]);

        let data = Data::try_from(vec![Field::new(
            Tag::ReadGroup,
            Value::String(String::from("rg0")),
        )])?;

        let header = Header::default();
        let record = Record::builder().set_data(data).build();

        writer.write_record(&header, &record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";
        assert_eq!(&writer.get_ref()[..], &expected[..]);

        Ok(())
    }
}
