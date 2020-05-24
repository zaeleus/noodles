mod alternative_allele;
mod builder;
mod contig;
mod filter;
pub mod format;
pub mod info;
mod number;
mod record;

pub use self::{
    alternative_allele::AlternativeAllele, builder::Builder, contig::Contig, filter::Filter,
    format::Format, info::Info, number::Number,
};

use std::{
    collections::HashMap,
    convert::TryFrom,
    error, fmt,
    str::{FromStr, Lines},
};

use self::record::Record;

static FILE_FORMAT: &str = "VCFv4.3";

#[derive(Debug)]
pub struct Header {
    file_format: String,
    infos: Vec<Info>,
    filters: Vec<Filter>,
    formats: Vec<Format>,
    alternative_alleles: Vec<AlternativeAllele>,
    assembly: Option<String>,
    contigs: Vec<Contig>,
    samples_names: Vec<String>,
    map: HashMap<String, String>,
}

impl Header {
    pub fn builder() -> Builder {
        Builder::new()
    }

    pub fn file_format(&self) -> &str {
        &self.file_format
    }

    pub fn infos(&self) -> &[Info] {
        &self.infos
    }

    pub fn filters(&self) -> &[Filter] {
        &self.filters
    }

    pub fn formats(&self) -> &[Format] {
        &self.formats
    }

    pub fn alternative_alleles(&self) -> &[AlternativeAllele] {
        &self.alternative_alleles
    }

    pub fn assembly(&self) -> Option<&str> {
        self.assembly.as_deref()
    }

    pub fn contigs(&self) -> &[Contig] {
        &self.contigs
    }

    pub fn sample_names(&self) -> &[String] {
        &self.samples_names
    }

    pub fn get(&self, key: &str) -> Option<&String> {
        self.map.get(key)
    }
}

impl Default for Header {
    fn default() -> Self {
        Builder::new().build()
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "##{}={}", record::Kind::FileFormat, self.file_format)?;

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
            writeln!(f, "##{}={}", record::Kind::Assembly, assembly)?;
        }

        for contig in self.contigs() {
            writeln!(f, "{}", contig)?;
        }

        for (key, value) in &self.map {
            writeln!(f, "##{}={}", key, value)?;
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

#[derive(Debug)]
pub enum ParseError {
    MissingFileFormat,
    InvalidRecord(record::ParseError),
    InvalidInfo(info::ParseError),
    InvalidFilter(filter::ParseError),
    InvalidFormat(format::ParseError),
    InvalidAlternativeAllele(alternative_allele::ParseError),
    InvalidContig(contig::ParseError),
    ExpectedEof,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("invalid header: ")?;

        match self {
            ParseError::MissingFileFormat => f.write_str("missing fileformat"),
            ParseError::InvalidRecord(e) => write!(f, "{}", e),
            ParseError::InvalidInfo(e) => write!(f, "{}", e),
            ParseError::InvalidFilter(e) => write!(f, "{}", e),
            ParseError::InvalidFormat(e) => write!(f, "{}", e),
            ParseError::InvalidAlternativeAllele(e) => write!(f, "{}", e),
            ParseError::InvalidContig(e) => write!(f, "{}", e),
            ParseError::ExpectedEof => f.write_str("expected EOF"),
        }
    }
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut header = Header::default();
        let mut lines = s.lines();

        header.file_format = parse_file_format(&mut lines)?;

        while let Some(line) = lines.next() {
            let record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record {
                Record::FileFormat(_) => todo!("unexpected fileformat"),
                Record::Info(fields) => {
                    let info = Info::try_from(&fields[..]).map_err(ParseError::InvalidInfo)?;
                    header.infos.push(info);
                }
                Record::Filter(fields) => {
                    let filter =
                        Filter::try_from(&fields[..]).map_err(ParseError::InvalidFilter)?;
                    header.filters.push(filter);
                }
                Record::Format(fields) => {
                    let format =
                        Format::try_from(&fields[..]).map_err(ParseError::InvalidFormat)?;
                    header.formats.push(format);
                }
                Record::AlternativeAllele(fields) => {
                    let alternative_allele = AlternativeAllele::try_from(&fields[..])
                        .map_err(ParseError::InvalidAlternativeAllele)?;
                    header.alternative_alleles.push(alternative_allele);
                }
                Record::Assembly(value) => {
                    header.assembly = Some(value);
                }
                Record::Contig(fields) => {
                    let contig =
                        Contig::try_from(&fields[..]).map_err(ParseError::InvalidContig)?;
                    header.contigs.push(contig);
                }
                Record::Header(samples_names) => {
                    header.samples_names = samples_names;
                    break;
                }
                Record::Other(key, value) => {
                    header.map.insert(key, value);
                }
            }
        }

        if lines.next().is_some() {
            return Err(ParseError::ExpectedEof);
        }

        Ok(header)
    }
}

fn parse_file_format(lines: &mut Lines) -> Result<String, ParseError> {
    let record = lines
        .next()
        .ok_or_else(|| ParseError::MissingFileFormat)
        .and_then(|line| line.parse().map_err(ParseError::InvalidRecord))?;

    match record {
        Record::FileFormat(value) => Ok(value),
        _ => Err(ParseError::MissingFileFormat),
    }
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
        let header = Header {
            file_format: FILE_FORMAT.into(),
            infos: vec![],
            filters: vec![],
            formats: vec![],
            alternative_alleles: vec![],
            assembly: Some(String::from("file:///assemblies.fasta")),
            contigs: vec![],
            samples_names: vec![],
            map: vec![(String::from("fileDate"), String::from("20200514"))]
                .into_iter()
                .collect(),
        };

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

        assert_eq!(header.get("fileDate"), Some(&String::from("20200506")));
        assert_eq!(header.get("source"), Some(&String::from("noodles-vcf")));

        Ok(())
    }

    #[test]
    fn test_from_str_without_file_format() {
        let s = r#"##ALT=<ID=DEL,Description="Deletion">
"#;

        assert!(matches!(
            s.parse::<Header>(),
            Err(ParseError::MissingFileFormat)
        ));
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

        assert!(matches!(s.parse::<Header>(), Err(ParseError::ExpectedEof)));
    }
}
