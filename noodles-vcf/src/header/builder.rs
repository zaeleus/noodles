use std::collections::HashMap;

use super::{AlternativeAllele, Contig, Filter, Format, Header, Info, FILE_FORMAT};

#[derive(Debug)]
pub struct Builder {
    file_format: String,
    infos: Vec<Info>,
    filters: Vec<Filter>,
    formats: Vec<Format>,
    alternative_alleles: Vec<AlternativeAllele>,
    assembly: Option<String>,
    contigs: Vec<Contig>,
    sample_names: Vec<String>,
    map: HashMap<String, String>,
}

impl Builder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_file_format<I>(mut self, file_format: I) -> Self
    where
        I: Into<String>,
    {
        self.file_format = file_format.into();
        self
    }

    pub fn add_info(mut self, info: Info) -> Self {
        self.infos.push(info);
        self
    }

    pub fn add_filter(mut self, filter: Filter) -> Self {
        self.filters.push(filter);
        self
    }

    pub fn add_format(mut self, format: Format) -> Self {
        self.formats.push(format);
        self
    }

    pub fn add_alternate_allele(mut self, alternative_allele: AlternativeAllele) -> Self {
        self.alternative_alleles.push(alternative_allele);
        self
    }

    pub fn set_assembly<I>(mut self, assembly: I) -> Self
    where
        I: Into<String>,
    {
        self.assembly = Some(assembly.into());
        self
    }

    pub fn add_contig(mut self, contig: Contig) -> Self {
        self.contigs.push(contig);
        self
    }

    pub fn add_sample_name<I>(mut self, sample_name: I) -> Self
    where
        I: Into<String>,
    {
        self.sample_names.push(sample_name.into());
        self
    }

    pub fn insert<K, V>(mut self, key: K, value: V) -> Self
    where
        K: Into<String>,
        V: Into<String>,
    {
        self.map.insert(key.into(), value.into());
        self
    }

    pub fn build(self) -> Header {
        Header {
            file_format: self.file_format,
            infos: self.infos,
            filters: self.filters,
            formats: self.formats,
            alternative_alleles: self.alternative_alleles,
            assembly: self.assembly,
            contigs: self.contigs,
            samples_names: self.sample_names,
            map: self.map,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            file_format: FILE_FORMAT.into(),
            infos: Vec::new(),
            filters: Vec::new(),
            formats: Vec::new(),
            alternative_alleles: Vec::new(),
            assembly: None,
            contigs: Vec::new(),
            sample_names: Vec::new(),
            map: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Builder::new().build();

        assert_eq!(header.file_format(), FILE_FORMAT);
        assert!(header.infos().is_empty());
        assert!(header.filters().is_empty());
        assert!(header.formats().is_empty());
        assert!(header.alternative_alleles().is_empty());
        assert!(header.assembly().is_none());
        assert!(header.contigs().is_empty());
        assert!(header.sample_names().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header::{format, info, Number},
            record::{self, alternate_bases::allele},
        };

        let header = Builder::new()
            .set_file_format("VCFv4.2")
            .add_info(Info::new(
                record::info::field::Key::SamplesWithDataCount,
                Number::Count(1),
                info::Type::Integer,
                String::from("Number of samples with data"),
            ))
            .add_filter(Filter::new(
                String::from("q10"),
                String::from("Quality below 10"),
            ))
            .add_format(Format::new(
                record::format::Key::Genotype,
                Number::Count(1),
                format::Type::String,
                String::from("Genotype"),
            ))
            .add_alternate_allele(AlternativeAllele::new(
                allele::Symbol::from(allele::symbol::Type::Deletion),
                String::from("Deletion"),
            ))
            .set_assembly("file:///assemblies.fasta")
            .add_contig(Contig::new(String::from("sq0")))
            .add_contig(Contig::new(String::from("sq1")))
            .add_sample_name("sample0")
            .insert("fileDate", "20200515")
            .build();

        assert_eq!(header.file_format(), "VCFv4.2");
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 2);
        assert_eq!(header.sample_names().len(), 1);
        assert_eq!(header.get("fileDate"), Some(&String::from("20200515")));

        Ok(())
    }
}
