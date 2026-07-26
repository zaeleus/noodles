mod line;
mod num;

use std::io::{self, Write};

use self::line::write_line;
use crate::{DirectiveBuf, LineBuf, feature::RecordBuf};

/// A GFF writer.
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W> {
    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(io::sink());
    /// let _inner = writer.get_ref();
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Returns a mutable reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    /// let mut writer = gff::io::Writer::new(io::sink());
    /// let _inner = writer.get_mut();
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(io::sink());
    /// let _inner = writer.into_inner();
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a GFF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gff as gff;
    /// let writer = gff::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Writes a [`LineBuf`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use bstr::BString;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}, LineBuf};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = LineBuf::Directive(gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// ));
    /// writer.write_line(&version)?;
    ///
    /// let comment = LineBuf::Comment(BString::from("noodles"));
    /// writer.write_line(&comment)?;
    ///
    /// let record = LineBuf::Record(gff::feature::RecordBuf::default());
    /// writer.write_line(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// #noodles
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    pub fn write_line(&mut self, line: &LineBuf) -> io::Result<()> {
        write_line(&mut self.inner, line)
    }

    /// Writes a GFF directive.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// assert_eq!(writer.get_ref(), b"##gff-version 3\n");
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_directive(&mut self, directive: &DirectiveBuf) -> io::Result<()> {
        line::write_directive(&mut self.inner, directive)?;
        line::write_newline(&mut self.inner)
    }

    /// Writes a GFF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// let record = gff::feature::RecordBuf::default();
    /// writer.write_record(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_record(&mut self, record: &RecordBuf) -> io::Result<()> {
        self.write_feature_record(record)
    }

    /// Writes a feature record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff::{self as gff, directive_buf::{key, Value}};
    ///
    /// let mut writer = gff::io::Writer::new(Vec::new());
    ///
    /// let version = gff::DirectiveBuf::new(
    ///     key::GFF_VERSION,
    ///     Some(Value::GffVersion(Default::default())),
    /// );
    /// writer.write_directive(&version)?;
    ///
    /// let record = gff::feature::RecordBuf::default();
    /// writer.write_feature_record(&record)?;
    ///
    /// let expected = b"##gff-version 3
    /// .\t.\t.\t1\t1\t.\t.\t.\t.
    /// ";
    ///
    /// assert_eq!(&writer.get_ref()[..], &expected[..]);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_feature_record(&mut self, record: &dyn crate::feature::Record) -> io::Result<()> {
        line::write_record(&mut self.inner, record)?;
        line::write_newline(&mut self.inner)
    }
}

#[cfg(test)]
mod tests {
    use bstr::BString;
    use noodles_core::Position;

    use super::*;
    use crate::io::Reader;

    // Values covering each character that must be encoded, plus ones that must not be.
    const VALUES: [&str; 10] = [
        "ndls",
        "ndls v1",
        "ndls\tv1",
        "ndls\nv1",
        "ndls\rv1",
        "ndls\x00v1",
        "ndls\x7fv1",
        "50%",
        "%25",
        "a;b=c&d,e",
    ];

    #[test]
    fn test_write_record_with_delimiters_in_fields() -> Result<(), Box<dyn std::error::Error>> {
        for value in VALUES {
            let record = RecordBuf::builder()
                .set_reference_sequence_name(BString::from(value))
                .set_source(BString::from(value))
                .set_type(BString::from(value))
                .set_start(Position::try_from(8)?)
                .set_end(Position::try_from(13)?)
                .build();

            let mut writer = Writer::new(Vec::new());
            writer.write_record(&record)?;
            let dst = writer.into_inner();

            // § "Description of the Format" (2020-08-18): "GFF3 files are nine-column,
            // tab-delimited, plain text files."
            let line_count = dst.iter().filter(|b| **b == b'\n').count();
            let field_count = dst.iter().filter(|b| **b == b'\t').count() + 1;
            assert_eq!((line_count, field_count), (1, 9), "{value:?}");

            let mut reader = Reader::new(&dst[..]);
            let actual = reader.record_bufs().next().transpose()?;
            assert_eq!(actual.as_ref(), Some(&record), "{value:?}");
        }

        Ok(())
    }

    #[test]
    fn test_write_feature_record_is_idempotent() -> Result<(), Box<dyn std::error::Error>> {
        for value in VALUES {
            let mut writer = Writer::new(Vec::new());

            let record = RecordBuf::builder()
                .set_reference_sequence_name(BString::from(value))
                .set_source(BString::from(value))
                .set_type(BString::from(value))
                .set_start(Position::try_from(8)?)
                .set_end(Position::try_from(13)?)
                .build();

            writer.write_record(&record)?;
            let expected = writer.into_inner();

            let mut reader = Reader::new(&expected[..]);
            let mut line = crate::Line::default();
            reader.read_line(&mut line)?;

            let record = line.as_record().transpose()?.expect("missing record");

            let mut writer = Writer::new(Vec::new());
            writer.write_feature_record(&record)?;

            assert_eq!(writer.into_inner(), expected, "{value:?}");
        }

        Ok(())
    }
}
