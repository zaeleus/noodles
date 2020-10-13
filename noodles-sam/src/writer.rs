use std::io::{self, Write};

use super::{record, Header, Record};

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
        let rnext = record
            .mate_reference_sequence_name()
            .map(|mate_reference_sequence_name| {
                if let Some(reference_sequence_name) = record.reference_sequence_name() {
                    if mate_reference_sequence_name == reference_sequence_name {
                        return "=";
                    }
                }

                mate_reference_sequence_name.as_str()
            })
            .unwrap_or(record::NULL_FIELD);

        let pos = record
            .position()
            .map(i32::from)
            .unwrap_or(record::position::UNMAPPED);

        let pnext = record
            .mate_position()
            .map(i32::from)
            .unwrap_or(record::position::UNMAPPED);

        write!(
            self.inner,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}",
            qname = record.read_name().as_ref(),
            flag = u16::from(record.flags()),
            rname = record.reference_sequence_name().map(|name| name.as_str()).unwrap_or(record::NULL_FIELD),
            pos = pos,
            mapq = u8::from(record.mapping_quality()),
            cigar = record.cigar(),
            rnext = rnext,
            pnext = pnext,
            tlen = record.template_len(),
            seq = record.sequence(),
            qual = record.quality_scores(),
        )?;

        if record.data().is_empty() {
            writeln!(self.inner)
        } else {
            writeln!(self.inner, "\t{}", record.data())
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::record::{
        data::{self, Field},
        Data,
    };

    use super::*;

    #[test]
    fn test_write_record_with_data() -> io::Result<()> {
        let mut writer = Writer::new(vec![]);

        let data = Data::from(vec![Field::new(
            data::field::Tag::ReadGroup,
            data::field::Value::String(String::from("rg0")),
        )]);

        let record = Record::builder().set_data(data).build();

        writer.write_record(&record)?;

        let expected = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";

        assert_eq!(&writer.get_ref()[..], &expected[..]);

        Ok(())
    }
}
