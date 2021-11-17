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
mod parser;
pub mod pedigree;
pub mod record;
pub mod sample;

pub use self::{
    alternative_allele::AlternativeAllele, builder::Builder, contig::Contig,
    file_format::FileFormat, filter::Filter, format::Format, info::Info, meta::Meta,
    number::Number, parser::ParseError, pedigree::Pedigree, record::Record, sample::Sample,
};

use std::str::FromStr;

use indexmap::{IndexMap, IndexSet};

/// VCF header info records.
pub type Infos = IndexMap<crate::record::info::field::Key, Info>;

/// VCF header filter records.
pub type Filters = IndexMap<String, Filter>;

/// VCF header format records.
pub type Formats = IndexMap<crate::record::genotypes::genotype::field::Key, Format>;

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
        self.file_format
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
    ///     .add_filter(Filter::new("q10", "Quality below 10"))
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
    /// use noodles_vcf::{self as vcf, header::Format, record::genotypes::genotype::field::Key};
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
    ///     .add_contig(Contig::new("sq0"))
    ///     .build();
    ///
    /// let contigs = header.contigs();
    /// assert_eq!(contigs.len(), 1);
    /// assert_eq!(contigs[0], Contig::new("sq0"));
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
    /// let expected: IndexSet<_> = [String::from("sample0"), String::from("sample1")]
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

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parser::parse(s)
    }
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
                [(String::from("Assay"), String::from("WholeGenome"))]
                    .into_iter()
                    .collect(),
            ))
            .add_pedigree(Pedigree::new(
                String::from("cid"),
                [
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
