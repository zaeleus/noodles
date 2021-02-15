use std::collections::HashMap;

use super::{AlternativeAllele, Contig, Filter, Format, Header, Info, Meta, Record, FILE_FORMAT};

/// A VCF header builder.
#[derive(Debug)]
pub struct Builder {
    file_format: String,
    infos: Vec<Info>,
    filters: Vec<Filter>,
    formats: Vec<Format>,
    alternative_alleles: Vec<AlternativeAllele>,
    assembly: Option<String>,
    contigs: Vec<Contig>,
    meta: Vec<Meta>,
    pedigree_db: Option<String>,
    sample_names: Vec<String>,
    map: HashMap<String, Record>,
}

impl Builder {
    /// Sets the fileformat record (`fileformat`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let header = vcf::Header::builder().set_file_format("VCFv4.3").build();
    /// assert_eq!(header.file_format(), "VCFv4.3");
    /// ```
    pub fn set_file_format<I>(mut self, file_format: I) -> Self
    where
        I: Into<String>,
    {
        self.file_format = file_format.into();
        self
    }

    /// Adds an information record (`INFO`).
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
    pub fn add_info(mut self, info: Info) -> Self {
        self.infos.push(info);
        self
    }

    /// Adds a filter record (`FILTER`).
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
    pub fn add_filter(mut self, filter: Filter) -> Self {
        self.filters.push(filter);
        self
    }

    /// Adds a genotype format record (`FORMAT`).
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
    pub fn add_format(mut self, format: Format) -> Self {
        self.formats.push(format);
        self
    }

    /// Adds an alternative allele record (`ALT`).
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
    pub fn add_alternative_allele(mut self, alternative_allele: AlternativeAllele) -> Self {
        self.alternative_alleles.push(alternative_allele);
        self
    }

    /// Sets an breakpoint assemblies record (`assembly`).
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
    pub fn set_assembly<I>(mut self, assembly: I) -> Self
    where
        I: Into<String>,
    {
        self.assembly = Some(assembly.into());
        self
    }

    /// Adds a contig record (`contig`).
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
    pub fn add_contig(mut self, contig: Contig) -> Self {
        self.contigs.push(contig);
        self
    }

    /// Adds a meta record (`META`).
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
    /// assert_eq!(header.meta(), [meta]);
    /// ```
    pub fn add_meta(mut self, meta: Meta) -> Self {
        self.meta.push(meta);
        self
    }

    /// Sets a pedigree database record (`pedigreeDB`).
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
    pub fn set_pedigree_db<I>(mut self, pedigree_db: I) -> Self
    where
        I: Into<String>,
    {
        self.pedigree_db = Some(pedigree_db.into());
        self
    }

    /// Adds a sample name.
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
    pub fn add_sample_name<I>(mut self, sample_name: I) -> Self
    where
        I: Into<String>,
    {
        self.sample_names.push(sample_name.into());
        self
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
    /// let header = vcf::Header::builder().insert(record.clone()).build();
    ///
    /// assert_eq!(header.get("fileDate"), Some(&record));
    /// ```
    pub fn insert(mut self, record: Record) -> Self {
        self.map.insert(record.key().to_string(), record);
        self
    }

    /// Builds a VCF header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf as vcf;
    /// let header = vcf::Header::builder().build();
    /// ```
    pub fn build(self) -> Header {
        Header {
            file_format: self.file_format,
            infos: self.infos,
            filters: self.filters,
            formats: self.formats,
            alternative_alleles: self.alternative_alleles,
            assembly: self.assembly,
            contigs: self.contigs,
            meta: self.meta,
            pedigree_db: self.pedigree_db,
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
            meta: Vec::new(),
            pedigree_db: None,
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
        let header = Builder::default().build();

        assert_eq!(header.file_format(), FILE_FORMAT);
        assert!(header.infos().is_empty());
        assert!(header.filters().is_empty());
        assert!(header.formats().is_empty());
        assert!(header.alternative_alleles().is_empty());
        assert!(header.assembly().is_none());
        assert!(header.contigs().is_empty());
        assert!(header.meta().is_empty());
        assert!(header.sample_names().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header::{self, format, info, Number},
            record::{self, alternate_bases::allele},
        };

        let record = Record::new(
            header::record::Key::Other(String::from("fileDate")),
            header::record::Value::String(String::from("20200709")),
        );

        let header = Builder::default()
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
                record::genotype::field::Key::Genotype,
                Number::Count(1),
                format::Type::String,
                String::from("Genotype"),
            ))
            .add_alternative_allele(AlternativeAllele::new(
                allele::Symbol::StructuralVariant(allele::symbol::StructuralVariant::from(
                    allele::symbol::structural_variant::Type::Deletion,
                )),
                String::from("Deletion"),
            ))
            .set_assembly("file:///assemblies.fasta")
            .add_contig(Contig::new(String::from("sq0")))
            .add_contig(Contig::new(String::from("sq1")))
            .add_meta(Meta::new(
                String::from("Assay"),
                vec![String::from("WholeGenome"), String::from("Exome")],
            ))
            .add_sample_name("sample0")
            .insert(record.clone())
            .build();

        assert_eq!(header.file_format(), "VCFv4.2");
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 2);
        assert_eq!(header.meta().len(), 1);
        assert_eq!(header.sample_names().len(), 1);
        assert_eq!(header.get("fileDate"), Some(&record));

        Ok(())
    }
}
