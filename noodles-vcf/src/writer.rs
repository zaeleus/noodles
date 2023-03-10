mod record;

use std::io::{self, Write};

use self::record::write_record;
use super::{Header, Record, VariantWriter};

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
/// let mut writer = vcf::Writer::new(Vec::new());
///
/// let header = vcf::Header::builder()
///     .add_contig("sq0".parse()?, Map::<Contig>::new())
///     .build();
///
/// writer.write_header(&header)?;
///
/// let record = vcf::Record::builder()
///     .set_chromosome("sq0".parse()?)
///     .set_position(Position::try_from(1)?)
///     .set_reference_bases("A".parse()?)
///     .build()?;
///
/// writer.write_record(&record);
///
/// let expected = b"##fileformat=VCFv4.3
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
    /// let writer = vcf::Writer::new(Vec::new());
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
    /// let writer = vcf::Writer::new(Vec::new());
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
    /// let mut writer = vcf::Writer::new(Vec::new());
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
    /// let writer = vcf::Writer::new(Vec::new());
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
    /// let mut writer = vcf::Writer::new(Vec::new());
    ///
    /// let header = vcf::Header::default();
    /// writer.write_header(&header)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn write_header(&mut self, header: &Header) -> io::Result<()> {
        write!(self.inner, "{header}")
    }

    /// Writes a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, record::Position};
    ///
    /// let record = vcf::Record::builder()
    ///     .set_chromosome("sq0".parse()?)
    ///     .set_position(Position::try_from(1)?)
    ///     .set_reference_bases("A".parse()?)
    ///     .build()?;
    ///
    /// let mut writer = vcf::Writer::new(Vec::new());
    /// writer.write_record(&record)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        write_record(&mut self.inner, record)
    }
}

impl<W> VariantWriter for Writer<W>
where
    W: Write,
{
    fn write_variant_header(&mut self, header: &Header) -> io::Result<()> {
        self.write_header(header)
    }

    fn write_variant_record(&mut self, _: &Header, record: &Record) -> io::Result<()> {
        self.write_record(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::record::Position;

    #[test]
    fn test_write_header() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let header = Header::default();
        writer.write_header(&header)?;

        let expected = b"##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(writer.get_ref().as_slice(), &expected[..]);

        Ok(())
    }

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .build()?;

        writer.write_record(&record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\n";

        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }

    #[test]
    fn test_write_record_with_format() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header::format::key,
            record::{
                genotypes::{genotype::field::Value, Keys},
                Genotypes,
            },
        };

        let mut writer = Writer::new(Vec::new());

        let genotypes = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![[
                (key::GENOTYPE, Some(Value::String(String::from("0|0")))),
                (key::CONDITIONAL_GENOTYPE_QUALITY, Some(Value::Integer(13))),
            ]
            .into_iter()
            .collect()],
        );

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_genotypes(genotypes)
            .build()?;

        writer.write_record(&record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13\n";

        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
