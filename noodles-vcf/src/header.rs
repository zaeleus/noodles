mod alternative_allele;
mod contig;
mod filter;
mod format;
mod info;
mod number;
mod record;

pub use self::{
    alternative_allele::AlternativeAllele, contig::Contig, filter::Filter, format::Format,
    info::Info, number::Number,
};

use std::{collections::HashMap, convert::TryFrom, str::FromStr};

use self::record::Record;

#[derive(Debug, Default)]
pub struct Header {
    file_format: String,
    infos: Vec<Info>,
    filters: Vec<Filter>,
    formats: Vec<Format>,
    alternative_alleles: Vec<AlternativeAllele>,
    assembly: Option<String>,
    contigs: Vec<Contig>,
    map: HashMap<String, String>,
}

impl Header {
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

    pub fn get(&self, key: &str) -> Option<&String> {
        self.map.get(key)
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
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut header = Header::default();
        let mut lines = s.lines();

        let record = lines
            .next()
            .ok_or_else(|| ParseError::MissingFileFormat)
            .and_then(|line| line.parse().map_err(ParseError::InvalidRecord))?;

        match record {
            Record::FileFormat(value) => header.file_format = value,
            _ => return Err(ParseError::MissingFileFormat),
        }

        for line in lines {
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
                Record::Other(key, value) => {
                    header.map.insert(key, value);
                }
            }
        }

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
"#;

        let header: Header = s.parse()?;

        assert_eq!(header.file_format(), "VCFv4.3");
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 3);

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
}
