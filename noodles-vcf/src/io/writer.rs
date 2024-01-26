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
/// use noodles_vcf::{
///     self as vcf,
///     header::record::value::{map::Contig, Map},
///     record::Position,
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
/// let record = vcf::Record::builder()
///     .set_chromosome("sq0")
///     .set_position(Position::from(1))
///     .set_reference_bases("A".parse()?)
///     .build()?;
///
/// writer.write_record(&header, &record);
///
/// let expected = b"##fileformat=VCFv4.4
/// ###contig=<ID=sq0>
/// #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
/// sq0\t1\t.\tA\t.\t.\t.\t.
/// ";
///
/// assert_eq!(&writer.get_ref()[..], &expected[..]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
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
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let header = vcf::Header::default();
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0")
    ///     .set_position(Position::from(1))
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let mut writer = vcf::io::Writer::new(Vec::new());
    /// writer.write_record(&header, &record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, _: &Header, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record)
    }
}

impl<W> crate::variant::io::Write for Writer<W>
where
    W: Write,
{
    fn write_variant_header(&mut self, header: &Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_variant_record(&mut self, header: &Header, record: &Record) -> io::Result<()> {
        self.write_record(header, record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Position;

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::default();

        let record = Record::builder()
            .set_chromosome("sq0")
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .build()?;

        let mut writer = Writer::new(Vec::new());
        writer.write_record(&header, &record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }

    #[test]
    fn test_write_record_with_format() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::{
            genotypes::{keys::key, sample::Value, Keys},
            Genotypes,
        };

        let header = Header::default();

        let genotypes = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![vec![
                Some(Value::String(String::from("0|0"))),
                Some(Value::Integer(13)),
            ]],
        );

        let record = Record::builder()
            .set_chromosome("sq0")
            .set_position(Position::from(1))
            .set_reference_bases("A".parse()?)
            .set_genotypes(genotypes)
            .build()?;

        let mut writer = Writer::new(Vec::new());
        writer.write_record(&header, &record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13\n";
        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
