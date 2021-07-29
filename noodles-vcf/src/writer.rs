use std::io::{self, Write};

use super::{Header, Record};

/// A VCF writer.
///
/// # Examples
///
/// ```
/// # use std::{convert::TryFrom, io};
/// use noodles_vcf::{self as vcf, header::Contig, record::Position};
///
/// let mut writer = vcf::Writer::new(Vec::new());
///
/// let header = vcf::Header::builder()
///     .add_contig(Contig::new("sq0"))
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
        write!(self.inner, "{}", header)
    }

    /// Writes a VCF record.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::convert::TryFrom;
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
        writeln!(self.inner, "{}", record)
    }
}

#[cfg(test)]
mod tests {
    use std::convert::TryFrom;

    use crate::record::{Format, Genotype, Position};

    use super::*;

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
        let mut writer = Writer::new(Vec::new());

        let format: Format = "GT:GQ".parse()?;

        let record = Record::builder()
            .set_chromosome("sq0".parse()?)
            .set_position(Position::try_from(1)?)
            .set_reference_bases("A".parse()?)
            .set_format(format.clone())
            .add_genotype(Genotype::from_str_format("0|0:13", &format)?)
            .build()?;

        writer.write_record(&record)?;

        let expected = b"sq0\t1\t.\tA\t.\t.\t.\t.\tGT:GQ\t0|0:13\n";

        assert_eq!(writer.get_ref(), expected);

        Ok(())
    }
}
