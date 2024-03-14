//! VCF writer.

mod builder;
mod header;
mod record;

use std::io::{self, Write};

pub use self::builder::Builder;
use self::{header::write_header, record::write_record};
use crate::{Header, Record};

/// A VCF writer.
///
/// # Examples
///
/// ```
/// # use std::io;
/// use noodles_core::Position;
/// use noodles_vcf::{
///     self as vcf,
///     header::record::value::{map::Contig, Map},
///     variant::io::Write,
/// };
///
/// let mut writer = vcf::io::Writer::new(Vec::new());
///
/// let header = vcf::Header::builder()
///     .add_contig("sq0", Map::<Contig>::new())
///     .build();
///
/// writer.write_header(&header)?;
///
/// let record = vcf::variant::RecordBuf::builder()
///     .set_reference_sequence_name("sq0")
///     .set_position(Position::MIN)
///     .set_reference_bases("A")
///     .build();
///
/// writer.write_variant_record(&header, &record);
///
/// let expected = b"##fileformat=VCFv4.4
/// ###contig=<ID=sq0>
/// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
/// sq0\t1\t.\tA\t.\t.\t.\t.
/// ";
///
/// assert_eq!(&writer.get_ref()[..], &expected[..]);
/// # Ok::<_, std::io::Error>(())
/// ```
#[derive(Debug)]
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: Write,
{
    /// Creates a VCF writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::io::Writer::new(Vec::new());
    /// ```
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::io::Writer::new(Vec::new());
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
    /// use noodles_vcf as vcf;
    /// let mut writer = vcf::io::Writer::new(Vec::new());
    /// assert!(writer.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Unwraps and returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let writer = vcf::io::Writer::new(Vec::new());
    /// assert!(writer.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Writes a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = vcf::io::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write_header(&mut self.inner, header)
    }

    /// Writes a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf as vcf;
    ///
    /// let mut writer = vcf::io::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// let record = vcf::Record::default();
    /// writer.write_record(&header, &record)?;
    ///
    /// assert_eq!(writer.get_ref(), b"sq0\t1\t.\tA\t.\t.\t.\t.\n");
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn write_record(&mut self, header: &Header, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, header, record)
    }
}

impl<W> crate::variant::io::Write for Writer<W>
where
    W: Write,
{
    fn write_variant_header(&mut self, header: &Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_variant_record(
        &mut self,
        header: &Header,
        record: &dyn crate::variant::Record,
    ) -> io::Result<()> {
        write_record(&mut self.inner, header, record)
    }
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::variant::{io::Write, RecordBuf};

    #[test]
    fn test_write_variant_record() -> io::Result<()> {
        let header = Header::default();

        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .build();

        let mut writer = Writer::new(Vec::new());
        writer.write_variant_record(&header, &record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }

    #[test]
    fn test_write_record_with_format() -> Result<(), Box<dyn std::error::Error>> {
        use crate::variant::{
            record::samples::keys::key,
            record_buf::{
                samples::{sample::Value, Keys},
                Samples,
            },
        };

        let header = Header::default();

        let samples = Samples::new(
            Keys::try_from(vec![
                String::from(key::GENOTYPE),
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
            ])?,
            vec![vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]],
        );

        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_position(Position::MIN)
            .set_reference_bases("A")
            .set_samples(samples)
            .build();

        let mut writer = Writer::new(Vec::new());
        writer.write_variant_record(&header, &record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
