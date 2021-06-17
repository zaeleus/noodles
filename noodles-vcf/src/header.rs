//! VCF header and fields.

pub mod alternative_allele;
mod builder;
pub mod contig;
pub mod file_format;
pub mod filter;
mod fmt;
pub mod format;
pub mod info;
pub mod meta;
mod number;
pub mod pedigree;
pub mod record;
pub mod sample;

pub use self::{
    alternative_allele::AlternativeAllele, builder::Builder, contig::Contig,
    file_format::FileFormat, filter::Filter, format::Format, info::Info, meta::Meta,
    number::Number, pedigree::Pedigree, record::Record, sample::Sample,
};

use std::{
    convert::TryFrom,
    error,
    str::{FromStr, Lines},
};

use indexmap::{IndexMap, IndexSet};

static HEADERS: &[&str] = &[
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
];
static FORMAT_HEADER: &str = "FORMAT";

/// VCF header info records.
pub type Infos = IndexMap<crate::record::info::field::Key, Info>;

/// VCF header filter records.
pub type Filters = IndexMap<String, Filter>;

/// VCF header format records.
pub type Formats = IndexMap<crate::record::genotype::field::Key, Format>;

/// VCF header alternative allele records.
pub type AlternativeAlleles =
    IndexMap<crate::record::alternate_bases::allele::Symbol, AlternativeAllele>;

/// VCF header contig records.
pub type Contigs = IndexMap<String, Contig>;

/// VCF header sample records.
pub type Samples = IndexMap<String, Sample>;

/// VCF header pedigree records.
pub type Pedigrees = IndexMap<String, Pedigree>;

/// VCF header sample names.
pub type SampleNames = IndexSet<String>;

/// A VCF header.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    file_format: FileFormat,
    infos: Infos,
    filters: Filters,
    formats: Formats,
    alternative_alleles: AlternativeAlleles,
    assembly: Option<String>,
    contigs: Contigs,
    meta: IndexMap<String, Meta>,
    samples: Samples,
    pedigrees: Pedigrees,
    pedigree_db: Option<String>,
    sample_names: SampleNames,
    map: IndexMap<String, Vec<Record>>,
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
    /// use noodles_vcf::{self as vcf, header::FileFormat};
    ///
    /// let header = vcf::Header::builder()
    ///     .set_file_format(FileFormat::default())
    ///     .build();
    ///
    /// assert_eq!(header.file_format(), FileFormat::default());
    /// ```
    pub fn file_format(&self) -> FileFormat {
        self.file_format.clone()
    }

    /// Returns a map of information records (`INFO`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Info, record::info::field::Key};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(Info::from(Key::SamplesWithDataCount))
    ///     .build();
    ///
    /// let infos = header.infos();
    /// assert_eq!(infos.len(), 1);
    /// assert_eq!(infos[0].id(), &Key::SamplesWithDataCount);
    /// ```
    pub fn infos(&self) -> &Infos {
        &self.infos
    }

    /// Returns a map of filter records (`FILTER`).
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
    pub fn filters(&self) -> &Filters {
        &self.filters
    }

    /// Returns a list of genotype format records (`FORMAT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Format, record::genotype::field::Key};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Format::from(Key::Genotype))
    ///     .build();
    ///
    /// let formats = header.formats();
    /// assert_eq!(formats.len(), 1);
    /// assert_eq!(formats[0].id(), &Key::Genotype);
    /// ```
    pub fn formats(&self) -> &Formats {
        &self.formats
    }

    /// Returns a map of symbolic alternate alleles (`ALT`).
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
    pub fn alternative_alleles(&self) -> &AlternativeAlleles {
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

    /// Returns a map of contig records (`contig`).
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
    /// let contigs = header.contigs();
    /// assert_eq!(contigs.len(), 1);
    /// assert_eq!(contigs[0], Contig::new(String::from("sq0")));
    /// ```
    pub fn contigs(&self) -> &Contigs {
        &self.contigs
    }

    /// Returns a map of meta records (`META`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Meta};
    ///
    /// let meta = Meta::new(
    ///     String::from("Assay"),
    ///     vec![String::from("WholeGenome"), String::from("Exome")],
    /// );
    ///
    /// let header = vcf::Header::builder()
    ///     .add_meta(meta.clone())
    ///     .build();
    ///
    /// let records = header.meta();
    /// assert_eq!(records.len(), 1);
    /// assert_eq!(records[0], meta);
    /// ```
    pub fn meta(&self) -> &IndexMap<String, Meta> {
        &self.meta
    }

    /// Returns a map of sample records (`SAMPLE`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Sample};
    ///
    /// let sample = Sample::new(String::from("sample0"), Default::default());
    ///
    /// let header = vcf::Header::builder()
    ///     .add_sample(sample.clone())
    ///     .build();
    ///
    /// let records = header.samples();
    /// assert_eq!(records.len(), 1);
    /// assert_eq!(records[0], sample);
    pub fn samples(&self) -> &Samples {
        &self.samples
    }

    /// Returns a map of pedigree records (`PEDIGREE`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Pedigree};
    ///
    /// let pedigree = Pedigree::new(
    ///     String::from("cid"),
    ///     vec![
    ///         (String::from("Father"), String::from("fid")),
    ///         (String::from("Mother"), String::from("mid")),
    ///     ]
    ///     .into_iter()
    ///     .collect(),
    /// );
    ///
    /// let header = vcf::Header::builder()
    ///     .add_pedigree(pedigree.clone())
    ///     .build();
    ///
    /// let records = header.pedigrees();
    /// assert_eq!(records.len(), 1);
    /// assert_eq!(records[0], pedigree);
    pub fn pedigrees(&self) -> &Pedigrees {
        &self.pedigrees
    }

    /// Returns a URI to the relationships between genomes (`pedigreeDB`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    ///
    /// let header = vcf::Header::builder()
    ///     .set_pedigree_db("file:///pedigree.db")
    ///     .build();
    ///
    /// assert_eq!(header.pedigree_db(), Some("file:///pedigree.db"));
    /// ```
    pub fn pedigree_db(&self) -> Option<&str> {
        self.pedigree_db.as_deref()
    }

    /// Returns a list sample names that come after the FORMAT column in the header record.
    ///
    /// # Examples
    ///
    /// ```
    /// use indexmap::IndexSet;
    /// use noodles_vcf as vcf;
    ///
    /// let header = vcf::Header::builder()
    ///     .add_sample_name("sample0")
    ///     .add_sample_name("sample1")
    ///     .build();
    ///
    /// let expected: IndexSet<_> = vec![String::from("sample0"), String::from("sample1")]
    ///     .into_iter()
    ///     .collect();
    ///
    /// assert_eq!(header.sample_names(), &expected);
    /// ```
    pub fn sample_names(&self) -> &SampleNames {
        &self.sample_names
    }

    /// Returns a header record with the given key.
    ///
    /// This includes all records other than `fileformat`, `INFO`, `FILTER`, `FORMAT`, `ALT`,
    /// `assembly`, `contig`, `META`, `SAMPLE`, and `pedigreeDB`.
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
    /// assert_eq!(header.get("fileDate"), Some(&[record][..]));
    /// assert_eq!(header.get("reference"), None);
    /// ```
    pub fn get(&self, key: &str) -> Option<&[Record]> {
        self.map.get(key).map(|r| &**r)
    }

    /// Inserts a key-value pair representing an unstructured record into the header.
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
    /// let mut header = vcf::Header::default();
    ///
    /// assert!(header.get("fileDate").is_none());
    ///
    /// header.insert(record.clone());
    ///
    /// assert_eq!(header.get("fileDate"), Some(&[record][..]));
    /// ```
    pub fn insert(&mut self, record: Record) {
        let key = record.key().to_string();
        let records = self.map.entry(key).or_default();
        records.push(record);
    }
}

impl Default for Header {
    fn default() -> Self {
        Builder::default().build()
    }
}

impl std::fmt::Display for Header {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}{}={}",
            record::PREFIX,
            record::Key::FileFormat,
            self.file_format()
        )?;

        for info in self.infos().values() {
            writeln!(f, "{}", info)?;
        }

        for filter in self.filters().values() {
            writeln!(f, "{}", filter)?;
        }

        for format in self.formats().values() {
            writeln!(f, "{}", format)?;
        }

        for alternative_allele in self.alternative_alleles().values() {
            writeln!(f, "{}", alternative_allele)?;
        }

        if let Some(assembly) = self.assembly() {
            writeln!(
                f,
                "{}{}={}",
                record::PREFIX,
                record::Key::Assembly,
                assembly
            )?;
        }

        for contig in self.contigs().values() {
            writeln!(f, "{}", contig)?;
        }

        for meta in self.meta().values() {
            writeln!(f, "{}", meta)?;
        }

        for sample in self.samples().values() {
            writeln!(f, "{}", sample)?;
        }

        for pedigree in self.pedigrees().values() {
            writeln!(f, "{}", pedigree)?;
        }

        if let Some(pedigree_db) = self.pedigree_db() {
            writeln!(
                f,
                "{}{}={}",
                record::PREFIX,
                record::Key::PedigreeDb,
                pedigree_db
            )?;
        }

        for records in self.map.values() {
            for record in records {
                writeln!(f, "{}{}={}", record::PREFIX, record.key(), record.value())?;
            }
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
    /// The file format (`fileformat`) is invalid.
    InvalidFileFormat(file_format::ParseError),
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
    /// A meta record (`META`) is invalid.
    InvalidMeta(meta::TryFromRecordError),
    /// A sample record (`SAMPLE`) is invalid.
    InvalidSample(sample::TryFromRecordError),
    /// A pedigree record (`PEDIGREE`) is invalid.
    InvalidPedigree(pedigree::TryFromRecordError),
    /// The header is missing.
    MissingHeader,
    /// The header is invalid.
    InvalidHeader(String, String),
    /// A sample name is duplicated.
    ///
    /// § 1.5 Header line syntax (2021-01-13): "Duplicate sample IDs are not allowed."
    DuplicateSampleName(String),
    /// More data unexpectedly appears after the header header (`#CHROM`...).
    ExpectedEof,
}

impl error::Error for ParseError {}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingFileFormat => f.write_str("missing fileformat"),
            Self::UnexpectedFileFormat => f.write_str("unexpected file format"),
            Self::InvalidFileFormat(e) => write!(f, "invalid file format: {}", e),
            Self::InvalidRecord(e) => write!(f, "invalid record: {}", e),
            Self::InvalidRecordValue => f.write_str("invalid record value"),
            Self::InvalidInfo(e) => write!(f, "invalid info: {}", e),
            Self::InvalidFilter(e) => write!(f, "invalid filter: {}", e),
            Self::InvalidFormat(e) => write!(f, "invalid format: {}", e),
            Self::InvalidAlternativeAllele(e) => {
                write!(f, "invalid alternative allele: {}", e)
            }
            Self::InvalidContig(e) => write!(f, "invalid contig: {}", e),
            Self::InvalidMeta(e) => write!(f, "invalid meta: {}", e),
            Self::InvalidSample(e) => write!(f, "invalid sample: {}", e),
            Self::InvalidPedigree(e) => write!(f, "invalid pedigree: {}", e),
            Self::MissingHeader => f.write_str("missing header"),
            Self::InvalidHeader(actual, expected) => {
                write!(f, "invalid header: expected {}, got {}", expected, actual)
            }
            Self::DuplicateSampleName(sample_name) => {
                write!(f, "duplicate sample name: {}", sample_name)
            }
            Self::ExpectedEof => f.write_str("expected EOF"),
        }
    }
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut builder = Self::builder();
        let mut lines = s.lines();

        let file_format = parse_file_format(&mut lines)?;
        builder = builder.set_file_format(file_format);

        let mut has_header = false;

        while let Some(line) = lines.next() {
            if line.starts_with("#CHROM") {
                builder = parse_header(builder, line)?;
                has_header = true;
                break;
            }

            builder = parse_record(builder, line)?;
        }

        if !has_header {
            return Err(ParseError::MissingHeader);
        }

        if lines.next().is_some() {
            return Err(ParseError::ExpectedEof);
        }

        Ok(builder.build())
    }
}

fn parse_file_format(lines: &mut Lines<'_>) -> Result<FileFormat, ParseError> {
    let record: Record = lines
        .next()
        .ok_or(ParseError::MissingFileFormat)
        .and_then(|line| line.parse().map_err(ParseError::InvalidRecord))?;

    if record.key() == &record::Key::FileFormat {
        match record.value() {
            record::Value::String(value) => value.parse().map_err(ParseError::InvalidFileFormat),
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
        record::Key::Meta => {
            let meta = Meta::try_from(record).map_err(ParseError::InvalidMeta)?;
            builder.add_meta(meta)
        }
        record::Key::Sample => {
            let sample = Sample::try_from(record).map_err(ParseError::InvalidSample)?;
            builder.add_sample(sample)
        }
        record::Key::Pedigree => {
            let pedigree = Pedigree::try_from(record).map_err(ParseError::InvalidPedigree)?;
            builder.add_pedigree(pedigree)
        }
        record::Key::PedigreeDb => match record.value() {
            record::Value::String(value) => builder.set_pedigree_db(value),
            _ => return Err(ParseError::InvalidRecordValue),
        },
        record::Key::Other(_) => builder.insert(record),
    };

    Ok(builder)
}

fn parse_header(mut builder: Builder, line: &str) -> Result<Builder, ParseError> {
    let mut fields = line.split(crate::record::FIELD_DELIMITER);

    for &expected in HEADERS.iter() {
        if let Some(actual) = fields.next() {
            if actual != expected {
                return Err(ParseError::InvalidHeader(actual.into(), expected.into()));
            }
        } else {
            return Err(ParseError::InvalidHeader(String::from(""), expected.into()));
        }
    }

    if let Some(field) = fields.next() {
        if field != FORMAT_HEADER {
            return Err(ParseError::InvalidHeader(
                field.into(),
                FORMAT_HEADER.into(),
            ));
        }

        let mut sample_names = IndexSet::new();

        for sample_name in fields {
            if !sample_names.insert(sample_name.into()) {
                return Err(ParseError::DuplicateSampleName(sample_name.into()));
            }
        }

        builder = builder.set_sample_names(sample_names);
    }

    Ok(builder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Header::default();
        assert_eq!(header.file_format(), FileFormat::default());
    }

    #[test]
    fn test_fmt() {
        let header = Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .add_filter(Filter::pass())
            .set_assembly("file:///assemblies.fasta")
            .add_meta(Meta::new(
                String::from("Assay"),
                vec![String::from("WholeGenome"), String::from("Exome")],
            ))
            .add_sample(Sample::new(
                String::from("sample0"),
                vec![(String::from("Assay"), String::from("WholeGenome"))]
                    .into_iter()
                    .collect(),
            ))
            .add_pedigree(Pedigree::new(
                String::from("cid"),
                vec![
                    (String::from("Father"), String::from("fid")),
                    (String::from("Mother"), String::from("mid")),
                ]
                .into_iter()
                .collect(),
            ))
            .insert(Record::new(
                record::Key::Other(String::from("fileDate")),
                record::Value::String(String::from("20200514")),
            ))
            .build();

        let expected = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##assembly=file:///assemblies.fasta
##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
##SAMPLE=<ID=sample0,Assay=WholeGenome>
##PEDIGREE=<ID=cid,Father=fid,Mother=mid>
##fileDate=20200514
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;

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
##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
##SAMPLE=<ID=sample0,Assay=WholeGenome>
##PEDIGREE=<ID=cid,Father=fid,Mother=mid>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let header: Header = s.parse()?;

        assert_eq!(header.file_format(), FileFormat::new(4, 3));
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 3);
        assert_eq!(header.meta().len(), 1);
        assert_eq!(header.samples().len(), 1);
        assert_eq!(header.pedigrees().len(), 1);
        assert_eq!(header.sample_names().len(), 1);

        assert_eq!(
            header.get("fileDate"),
            Some(
                &[Record::new(
                    record::Key::Other(String::from("fileDate")),
                    record::Value::String(String::from("20200506")),
                )][..]
            )
        );

        assert_eq!(
            header.get("source"),
            Some(
                &[Record::new(
                    record::Key::Other(String::from("source")),
                    record::Value::String(String::from("noodles-vcf")),
                )][..]
            )
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
        let s = r#"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;
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

    #[test]
    fn test_from_str_with_missing_headers() {
        let s = "##fileformat=VCFv4.3
";
        assert_eq!(s.parse::<Header>(), Err(ParseError::MissingHeader));
    }

    #[test]
    fn test_from_str_with_invalid_headers() {
        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUALITY	FILTER	INFO
";

        assert_eq!(
            s.parse::<Header>(),
            Err(ParseError::InvalidHeader(
                String::from("QUALITY"),
                String::from("QUAL")
            ))
        );

        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID
";

        assert_eq!(
            s.parse::<Header>(),
            Err(ParseError::InvalidHeader(
                String::from(""),
                String::from("REF")
            ))
        );

        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	sample0
";

        assert_eq!(
            s.parse::<Header>(),
            Err(ParseError::InvalidHeader(
                String::from("sample0"),
                String::from("FORMAT")
            ))
        );
    }

    #[test]
    fn test_from_str_with_duplicate_sample_names() {
        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0	sample0
";

        assert_eq!(
            s.parse::<Header>(),
            Err(ParseError::DuplicateSampleName(String::from("sample0")))
        );
    }

    #[test]
    fn test_insert_with_duplicate_keys() {
        let records = [
            Record::new(
                record::Key::Other(String::from("noodles")),
                record::Value::String(String::from("0")),
            ),
            Record::new(
                record::Key::Other(String::from("noodles")),
                record::Value::String(String::from("1")),
            ),
        ];

        let mut header = Header::default();

        for record in &records {
            header.insert(record.clone());
        }

        assert_eq!(header.get("noodles"), Some(&records[..]));
    }
}
