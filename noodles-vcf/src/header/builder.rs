use super::{
    record, AlternativeAllele, AlternativeAlleles, Contig, Contigs, FileFormat, Filter, Filters,
    Format, Formats, Header, Info, Infos, Meta, Pedigree, Pedigrees, Records, Sample, SampleNames,
    Samples,
};

use indexmap::IndexMap;

/// A VCF header builder.
#[derive(Debug, Default)]
pub struct Builder {
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
    other_records: Records,
}

impl Builder {
    /// Sets the fileformat record (`fileformat`).
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
    pub fn set_file_format(mut self, file_format: FileFormat) -> Self {
        self.file_format = file_format;
        self
    }

    /// Adds an information record (`INFO`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::{info::Key, Info}};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(Info::from(Key::SamplesWithDataCount))
    ///     .build();
    ///
    /// let infos = header.infos();
    /// assert_eq!(infos.len(), 1);
    /// assert_eq!(infos[0].id(), &Key::SamplesWithDataCount);
    /// ```
    pub fn add_info(mut self, info: Info) -> Self {
        self.infos.insert(info.id().clone(), info);
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
    ///     .add_filter(Filter::new("q10", "Quality below 10"))
    ///     .build();
    ///
    /// let filters = header.filters();
    /// assert_eq!(filters.len(), 1);
    /// assert_eq!(filters[0].id(), "q10");
    /// ```
    pub fn add_filter(mut self, filter: Filter) -> Self {
        self.filters.insert(filter.id().into(), filter);
        self
    }

    /// Adds a genotype format record (`FORMAT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::{format::Key, Format}};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(Format::from(Key::Genotype))
    ///     .build();
    ///
    /// let formats = header.formats();
    /// assert_eq!(formats.len(), 1);
    /// assert_eq!(formats[0].id(), &Key::Genotype);
    /// ```
    pub fn add_format(mut self, format: Format) -> Self {
        self.formats.insert(format.id().clone(), format);
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
    ///         "Deletion",
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
        self.alternative_alleles
            .insert(alternative_allele.id().clone(), alternative_allele);
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
    ///     .add_contig(Contig::new("sq0".parse()?))
    ///     .build();
    ///
    /// let contigs = header.contigs();
    /// assert_eq!(contigs.len(), 1);
    /// assert_eq!(contigs[0], Contig::new("sq0".parse()?));
    /// # Ok::<_, vcf::header::contig::name::ParseError>(())
    /// ```
    pub fn add_contig(mut self, contig: Contig) -> Self {
        self.contigs.insert(contig.id().as_ref().into(), contig);
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
    /// let records = header.meta();
    /// assert_eq!(records.len(), 1);
    /// assert_eq!(records[0], meta);
    /// ```
    pub fn add_meta(mut self, meta: Meta) -> Self {
        self.meta.insert(meta.id().into(), meta);
        self
    }

    /// Adds a sample record (`SAMPLE`).
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
    /// ```
    pub fn add_sample(mut self, sample: Sample) -> Self {
        self.samples.insert(sample.id().into(), sample);
        self
    }

    /// Adds a pedigree record (`PEDIGREE`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::Pedigree};
    ///
    /// let pedigree = Pedigree::new(
    ///     String::from("cid"),
    ///     [
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
    /// ```
    pub fn add_pedigree(mut self, pedigree: Pedigree) -> Self {
        self.pedigrees.insert(pedigree.id().into(), pedigree);
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

    /// Sets sample names.
    ///
    /// # Examples
    ///
    /// ```
    /// use indexmap::IndexSet;
    /// use noodles_vcf as vcf;
    ///
    /// let sample_names: IndexSet<_> = [String::from("sample0")]
    ///     .into_iter()
    ///     .collect();
    ///
    /// let header = vcf::Header::builder()
    ///     .set_sample_names(sample_names.clone())
    ///     .build();
    ///
    /// assert_eq!(header.sample_names(), &sample_names);
    /// ```
    pub fn set_sample_names(mut self, sample_names: SampleNames) -> Self {
        self.sample_names = sample_names;
        self
    }

    /// Adds a sample name.
    ///
    /// Duplicate names are discarded.
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
    /// let expected: IndexSet<_> = [String::from("sample0"), String::from("sample1")]
    ///     .into_iter()
    ///     .collect();
    ///
    /// assert_eq!(header.sample_names(), &expected);
    /// ```
    pub fn add_sample_name<I>(mut self, sample_name: I) -> Self
    where
        I: Into<String>,
    {
        self.sample_names.insert(sample_name.into());
        self
    }

    /// Inserts a key-value pair representing an unstructured record into the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::{record::{Key, Value}}};
    /// let (key, value) = (Key::from("fileDate"), Value::from("20200709"));
    /// let header = vcf::Header::builder().insert(key.clone(), value.clone()).build();
    /// assert_eq!(header.get(&key), Some(&[value][..]));
    /// ```
    pub fn insert(mut self, key: record::Key, value: record::Value) -> Self {
        let records = self.other_records.entry(key).or_default();
        records.push(value);
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
            samples: self.samples,
            pedigrees: self.pedigrees,
            pedigree_db: self.pedigree_db,
            sample_names: self.sample_names,
            other_records: self.other_records,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let header = Builder::default().build();

        assert_eq!(header.file_format(), FileFormat::default());
        assert!(header.infos().is_empty());
        assert!(header.filters().is_empty());
        assert!(header.formats().is_empty());
        assert!(header.alternative_alleles().is_empty());
        assert!(header.assembly().is_none());
        assert!(header.contigs().is_empty());
        assert!(header.meta().is_empty());
        assert!(header.samples().is_empty());
        assert!(header.pedigrees().is_empty());
        assert!(header.pedigree_db().is_none());
        assert!(header.sample_names().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), crate::header::contig::name::ParseError> {
        use crate::{
            header::{self, format::Key as FormatKey, info::Key as InfoKey},
            record::alternate_bases::allele,
        };

        let (key, value) = (
            header::record::Key::from("fileDate"),
            header::record::Value::from("20200709"),
        );

        let header = Builder::default()
            .set_file_format(FileFormat::new(4, 3))
            .add_info(Info::from(InfoKey::SamplesWithDataCount))
            .add_filter(Filter::new("q10", "Quality below 10"))
            .add_format(Format::from(FormatKey::Genotype))
            .add_alternative_allele(AlternativeAllele::new(
                allele::Symbol::StructuralVariant(allele::symbol::StructuralVariant::from(
                    allele::symbol::structural_variant::Type::Deletion,
                )),
                "Deletion",
            ))
            .set_assembly("file:///assemblies.fasta")
            .add_contig(Contig::new("sq0".parse()?))
            .add_contig(Contig::new("sq1".parse()?))
            .add_meta(Meta::new(
                String::from("Assay"),
                vec![String::from("WholeGenome"), String::from("Exome")],
            ))
            .add_sample(Sample::new(String::from("sample0"), Default::default()))
            .add_sample_name("sample0")
            .add_pedigree(Pedigree::new(
                String::from("cid"),
                [
                    (String::from("Father"), String::from("fid")),
                    (String::from("Mother"), String::from("mid")),
                ]
                .into_iter()
                .collect(),
            ))
            .insert(key.clone(), value.clone())
            .insert(key.clone(), value.clone())
            .build();

        assert_eq!(header.file_format(), FileFormat::new(4, 3));
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 2);
        assert_eq!(header.meta().len(), 1);
        assert_eq!(header.samples().len(), 1);
        assert_eq!(header.sample_names().len(), 1);
        assert_eq!(header.get(&key), Some(&[value.clone(), value][..]));

        Ok(())
    }
}
