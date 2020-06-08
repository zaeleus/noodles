use std::io::{self, Write};

use super::{Header, Record};

/// A SAM writer.
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
    /// assert_eq!(writer.get_ref(), b"*\t0\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write!(
            self.inner,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}",
            qname = record.read_name().as_ref(),
            flag = u16::from(record.flags()),
            rname = record.reference_sequence_name().as_ref(),
            pos = i32::from(record.position()),
            mapq = u8::from(record.mapping_quality()),
            cigar = record.cigar(),
            rnext = record.mate_reference_sequence_name().as_ref(),
            pnext = i32::from(record.mate_position()),
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
    use crate::{
        header,
        record::{
            data::{self, Field},
            Data,
        },
    };

    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Writer::new(vec![]);

        let header = Header::builder()
            .set_header(header::header::Header::default())
            .build();

        writer.write_header(&header)?;

        let expected = b"@HD\tVN:1.6\n";

        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut writer = Writer::new(vec![]);

        let record = Record::default();
        writer.write_record(&record)?;

        let expected = b"*\t0\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";

        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }

    #[test]
    fn test_write_record_with_data() -> io::Result<()> {
        let mut writer = Writer::new(vec![]);

        let data = Data::new(vec![Field::new(
            data::field::Tag::ReadGroup,
            data::field::Value::String(String::from("rg0")),
        )]);

        let record = Record::builder().set_data(data).build();

        writer.write_record(&record)?;

        let expected = b"*\t0\t*\t0\t255\t*\t*\t0\t0\t*\t*\tRG:Z:rg0\n";

        assert_eq!(&writer.get_ref()[..], &expected[..]);

        Ok(())
    }
}
