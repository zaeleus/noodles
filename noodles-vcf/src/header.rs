//! VCF header and fields.

mod alternative_allele;
mod builder;
mod contig;
mod filter;
pub mod format;
pub mod info;
mod number;
pub mod record;

pub use self::{
    alternative_allele::AlternativeAllele, builder::Builder, contig::Contig, filter::Filter,
    format::Format, info::Info, number::Number, record::Record,
};

use std::{
    collections::HashMap,
    convert::TryFrom,
    error, fmt,
    str::{FromStr, Lines},
};

static FILE_FORMAT: &str = "VCFv4.3";

/// A VCF header.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    file_format: String,
    infos: Vec<Info>,
    filters: Vec<Filter>,
    formats: Vec<Format>,
    alternative_alleles: Vec<AlternativeAllele>,
    assembly: Option<String>,
    contigs: Vec<Contig>,
    samples_names: Vec<String>,
    map: HashMap<String, Record>,
}

impl Header {
    /// Returns a builder to create a record from each of its fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let builder = vcf::Header::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the file format (`fileformat`) of the VCF.
    ///
    /// `fileformat` is a reqiured meta record and is guaranteed to be set.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let header = vcf::Header::builder().set_file_format("VCFv4.3").build();
    /// assert_eq!(header.file_format(), "VCFv4.3");
    /// ```
    pub fn file_format(&self) -> &str {
        &self.file_format
    }

    /// Returns a list of information records (`INFO`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{info::Type, Info, Number},
    ///     record::info::field::Key,
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(Info::new(
    ///         Key::SamplesWithDataCount,
    ///         Number::Count(1),
    ///         Type::Integer,
    ///         String::from("Number of samples with data"),
    ///     ))
    ///     .build();
    ///
    /// let infos = header.infos();
    /// assert_eq!(infos.len(), 1);
    /// assert_eq!(infos[0].id(), &Key::SamplesWithDataCount);
    /// ```
    pub fn infos(&self) -> &[Info] {
        &self.infos
    }

    /// Returns a list of filter records (`FILTER`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Filter};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_filter(Filter::new(
    ///         String::from("q10"),
    ///         String::from("Quality below 10"),
    ///     ))
    ///     .build();
    ///
    /// let filters = header.filters();
    /// assert_eq!(filters.len(), 1);
    /// assert_eq!(filters[0].id(), "q10");
    /// ```
    pub fn filters(&self) -> &[Filter] {
        &self.filters
    }

    /// Returns a list of genotype format records (`FORMAT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Type, Format, Number},
    ///     record::genotype::field::Key,
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Format::new(
    ///         Key::Genotype,
    ///         Number::Count(1),
    ///         Type::String,
    ///         String::from("Genotype"),
    ///     ))
    ///     .build();
    ///
    /// let formats = header.formats();
    /// assert_eq!(formats.len(), 1);
    /// assert_eq!(formats[0].id(), &Key::Genotype);
    /// ```
    pub fn formats(&self) -> &[Format] {
        &self.formats
    }

    /// Returns a list of symbolic alternate alleles (`ALT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::AlternativeAllele,
    ///     record::alternate_bases::allele::{
    ///         symbol::{structural_variant::Type, StructuralVariant},
    ///         Symbol,
    ///     },
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_alternative_allele(AlternativeAllele::new(
    ///         Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
    ///         String::from("Deletion"),
    ///     ))
    ///     .build();
    ///
    /// let alternative_alleles = header.alternative_alleles();
    /// assert_eq!(alternative_alleles.len(), 1);
    /// assert_eq!(
    ///     alternative_alleles[0].id(),
    ///     &Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion))
    /// );
    /// ```
    pub fn alternative_alleles(&self) -> &[AlternativeAllele] {
        &self.alternative_alleles
    }

    /// Returns a URI to the breakpoint assemblies (`assembly`) referenced in records.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let header = vcf::Header::builder()
    ///     .set_assembly("file:///assemblies.fasta")
    ///     .build();
    ///
    /// assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
    /// ```
    pub fn assembly(&self) -> Option<&str> {
        self.assembly.as_deref()
    }

    /// Returns a list of contig records (`contig`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Contig};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_contig(Contig::new(String::from("sq0")))
    ///     .build();
    ///
    /// assert_eq!(header.contigs(), [Contig::new(String::from("sq0"))]);
    /// ```
    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    /// Returns a list sample names that come after the FORMAT column in the header record.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let header = vcf::Header::builder()
    ///     .add_sample_name("sample0")
    ///     .add_sample_name("sample1")
    ///     .build();
    ///
    /// assert_eq!(header.sample_names(), [
    ///     String::from("sample0"),
    ///     String::from("sample1"),
    /// ]);
    /// ```
    pub fn sample_names(&self) -> &[String] {
        &self.samples_names
    }

    /// Returns a header record with the given key.
    ///
    /// This includes all records other than `fileformat`, `INFO`, `FILTER`, `FORMAT`, `ALT`,
    /// `assembly`, and `contig`.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::{record::{Key, Value}, Record}};
    ///
    /// let record = Record::new(
    ///     Key::Other(String::from("fileDate")),
    ///     Value::String(String::from("20200709")),
    /// );
    ///
    /// let header = vcf::Header::builder().insert(record.clone()).build();
    ///
    /// assert_eq!(header.get("fileDate"), Some(&record));
    /// assert_eq!(header.get("reference"), None);
    /// ```
    pub fn get(&self, key: &str) -> Option<&Record> {
        self.map.get(key)
    }
}

impl Default for Header {
    fn default() -> Self {
        Builder::default().build()
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "##{}={}", record::Key::FileFormat, self.file_format)?;

        for info in self.infos() {
            writeln!(f, "{}", info)?;
        }

        for format in self.formats() {
            writeln!(f, "{}", format)?;
        }

        for alternative_allele in self.alternative_alleles() {
            writeln!(f, "{}", alternative_allele)?;
        }

        if let Some(assembly) = self.assembly() {
            writeln!(f, "##{}={}", record::Key::Assembly, assembly)?;
        }

        for contig in self.contigs() {
            writeln!(f, "{}", contig)?;
        }

        for record in self.map.values() {
            writeln!(f, "##{}={}", record.key(), record.value())?;
        }

        f.write_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

        if !self.sample_names().is_empty() {
            f.write_str("\tFORMAT")?;

            for sample_name in self.sample_names() {
                write!(f, "\t{}", sample_name)?;
            }
        }

        f.write_str("\n")?;

        Ok(())
    }
}

/// An error returned when a raw VCF header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The file format (`fileformat`) is missing.
    MissingFileFormat,
    /// The file format (`fileformat`) appears other than the first line.
    UnexpectedFileFormat,
    /// A record is invalid.
    InvalidRecord(record::ParseError),
    /// A record has an invalid value.
    InvalidRecordValue,
    /// An information record (`INFO`) is invalid.
    InvalidInfo(info::TryFromRecordError),
    /// A filter record (`FILTER`) is invalid.
    InvalidFilter(filter::TryFromRecordError),
    /// A genotype format record (`FORMAT`) is invalid.
    InvalidFormat(format::TryFromRecordError),
    /// A symboloic alternate allele record (`ALT`) is invalid.
    InvalidAlternativeAllele(alternative_allele::TryFromRecordError),
    /// A contig record (`contig`) is invalid.
    InvalidContig(contig::TryFromRecordError),
    /// More data unexpectedly appears after the header header (`#CHROM`...).
    ExpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::MissingFileFormat => f.write_str("missing fileformat"),
            ParseError::UnexpectedFileFormat => f.write_str("unexpected file format"),
            ParseError::InvalidRecord(e) => write!(f, "invalid record: {}", e),
            ParseError::InvalidRecordValue => f.write_str("invalid record value"),
            ParseError::InvalidInfo(e) => write!(f, "invalid info: {}", e),
            ParseError::InvalidFilter(e) => write!(f, "invalid filter: {}", e),
            ParseError::InvalidFormat(e) => write!(f, "invalid format: {}", e),
            ParseError::InvalidAlternativeAllele(e) => {
                write!(f, "invalid alternative allele: {}", e)
            }
            ParseError::InvalidContig(e) => write!(f, "invalid contig: {}", e),
            ParseError::ExpectedEof => f.write_str("expected EOF"),
        }
    }
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut builder = Header::builder();
        let mut lines = s.lines();

        let file_format = parse_file_format(&mut lines)?;
        builder = builder.set_file_format(file_format);

        while let Some(line) = lines.next() {
            if line.starts_with("#CHROM") {
                builder = parse_header(builder, line)?;
                break;
            } else {
                builder = parse_record(builder, line)?;
            }
        }

        if lines.next().is_some() {
            return Err(ParseError::ExpectedEof);
        }

        Ok(builder.build())
    }
}

fn parse_file_format(lines: &mut Lines<'_>) -> Result<String, ParseError> {
    let record: Record = lines
        .next()
        .ok_or_else(|| ParseError::MissingFileFormat)
        .and_then(|line| line.parse().map_err(ParseError::InvalidRecord))?;

    if record.key() == &record::Key::FileFormat {
        match record.value() {
            record::Value::String(value) => Ok(value.into()),
            _ => Err(ParseError::InvalidRecordValue),
        }
    } else {
        Err(ParseError::MissingFileFormat)
    }
}

fn parse_record(mut builder: Builder, line: &str) -> Result<Builder, ParseError> {
    let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;

    builder = match record.key() {
        record::Key::FileFormat => {
            return Err(ParseError::UnexpectedFileFormat);
        }
        record::Key::Info => {
            let info = Info::try_from(record).map_err(ParseError::InvalidInfo)?;
            builder.add_info(info)
        }
        record::Key::Filter => {
            let filter = Filter::try_from(record).map_err(ParseError::InvalidFilter)?;
            builder.add_filter(filter)
        }
        record::Key::Format => {
            let format = Format::try_from(record).map_err(ParseError::InvalidFormat)?;
            builder.add_format(format)
        }
        record::Key::AlternativeAllele => {
            let alternative_allele = AlternativeAllele::try_from(record)
                .map_err(ParseError::InvalidAlternativeAllele)?;
            builder.add_alternative_allele(alternative_allele)
        }
        record::Key::Assembly => match record.value() {
            record::Value::String(value) => builder.set_assembly(value),
            _ => return Err(ParseError::InvalidRecordValue),
        },
        record::Key::Contig => {
            let contig = Contig::try_from(record).map_err(ParseError::InvalidContig)?;
            builder.add_contig(contig)
        }
        record::Key::Other(_) => builder.insert(record),
    };

    Ok(builder)
}

fn parse_header(mut builder: Builder, line: &str) -> Result<Builder, ParseError> {
    let sample_names = line.split(crate::record::FIELD_DELIMITER).skip(9);

    for sample_name in sample_names {
        builder = builder.add_sample_name(sample_name);
    }

    Ok(builder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Header::default();
        assert_eq!(header.file_format(), FILE_FORMAT);
    }

    #[test]
    fn test_fmt() {
        let header = Header::builder()
            .set_file_format("VCFv4.3")
            .set_assembly("file:///assemblies.fasta")
            .insert(Record::new(
                record::Key::Other(String::from("fileDate")),
                record::Value::String(String::from("20200514")),
            ))
            .build();

        let expected = "\
##fileformat=VCFv4.3
##assembly=file:///assemblies.fasta
##fileDate=20200514
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

        assert_eq!(header.to_string(), expected);
    }

    #[test]
    fn test_fmt_with_genotypes() {
        let header = Header::builder().add_sample_name("sample0").build();
        let expected = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0
";
        assert_eq!(header.to_string(), expected);

        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();

        let expected = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample0\tsample1
";

        assert_eq!(header.to_string(), expected);
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20200506
##source=noodles-vcf
##assembly=file:///assemblies.fasta
##contig=<ID=sq0,length=8>
##contig=<ID=sq1,length=13>
##contig=<ID=sq2,length=21>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##FILTER=<ID=q10,Description="Quality below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;

        let header: Header = s.parse()?;

        assert_eq!(header.file_format(), "VCFv4.3");
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 3);
        assert!(header.sample_names().is_empty());

        assert_eq!(
            header.get("fileDate"),
            Some(&Record::new(
                record::Key::Other(String::from("fileDate")),
                record::Value::String(String::from("20200506")),
            ))
        );

        assert_eq!(
            header.get("source"),
            Some(&Record::new(
                record::Key::Other(String::from("source")),
                record::Value::String(String::from("noodles-vcf")),
            ))
        );

        Ok(())
    }

    #[test]
    fn test_from_str_without_file_format() {
        let s = r#"##ALT=<ID=DEL,Description="Deletion">
"#;

        assert_eq!(s.parse::<Header>(), Err(ParseError::MissingFileFormat));
    }

    #[test]
    fn test_from_str_without_assembly() -> Result<(), ParseError> {
        let s = r#"##fileformat=VCFv4.3"#;
        let header: Header = s.parse()?;
        assert!(header.assembly().is_none());
        Ok(())
    }

    #[test]
    fn test_from_str_with_data_after_header() {
        let s = r#"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
##contig=<ID=sq0,length=8>
"#;

        assert_eq!(s.parse::<Header>(), Err(ParseError::ExpectedEof));
    }

    #[test]
    fn test_from_str_with_multiple_fileformats() {
        let s = "\
##fileformat=VCFv4.3
##fileformat=VCFv4.3
";

        assert_eq!(s.parse::<Header>(), Err(ParseError::UnexpectedFileFormat));
    }
}
