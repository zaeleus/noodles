use super::{
    record::{
        self,
        value::{
            map::{AlternativeAllele, Contig, Filter, Format, Info},
            Map,
        },
    },
    AlternativeAlleles, Contigs, FileFormat, Filters, Formats, Header, Infos, OtherRecords,
    SampleNames, StringMaps,
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
    contigs: Contigs,
    sample_names: SampleNames,
    other_records: OtherRecords,
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
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map::Info, Map},
    ///     variant::record::info::field::key,
    /// };
    ///
    /// let id = key::SAMPLES_WITH_DATA_COUNT;
    /// let info = Map::<Info>::from(id);
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(id, info.clone())
    ///     .build();
    ///
    /// let infos = header.infos();
    /// assert_eq!(infos.len(), 1);
    /// assert_eq!(&infos[0], &info);
    /// ```
    pub fn add_info<I>(mut self, id: I, info: Map<Info>) -> Self
    where
        I: Into<String>,
    {
        self.infos.insert(id.into(), info);
        self
    }

    /// Adds a filter record (`FILTER`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::record::value::{map::Filter, Map}};
    ///
    /// let filter = Map::<Filter>::new("Quality below 10");
    ///
    /// let header = vcf::Header::builder()
    ///     .add_filter("q10", filter.clone())
    ///     .build();
    ///
    /// let filters = header.filters();
    /// assert_eq!(filters.len(), 1);
    /// assert_eq!(&filters[0], &filter);
    /// ```
    pub fn add_filter<I>(mut self, id: I, filter: Map<Filter>) -> Self
    where
        I: Into<String>,
    {
        self.filters.insert(id.into(), filter);
        self
    }

    /// Adds a genotype format record (`FORMAT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map::Format, Map},
    ///     variant::record_buf::samples::keys::key,
    /// };
    ///
    /// let id = key::GENOTYPE;
    /// let format = Map::<Format>::from(id);
    ///
    /// let header = vcf::Header::builder()
    ///     .add_format(id, format.clone())
    ///     .build();
    ///
    /// let formats = header.formats();
    /// assert_eq!(formats.len(), 1);
    /// assert_eq!(&formats[0], &format);
    /// ```
    pub fn add_format<I>(mut self, id: I, format: Map<Format>) -> Self
    where
        I: Into<String>,
    {
        self.formats.insert(id.into(), format);
        self
    }

    /// Adds an alternative allele record (`ALT`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map::AlternativeAllele, Map},
    /// };
    ///
    /// let alt = Map::<AlternativeAllele>::new("Deletion");
    ///
    /// let header = vcf::Header::builder()
    ///     .add_alternative_allele("DEL", alt.clone())
    ///     .build();
    ///
    /// let alternative_alleles = header.alternative_alleles();
    /// assert_eq!(alternative_alleles.len(), 1);
    /// assert_eq!(&alternative_alleles[0], &alt);
    /// ```
    pub fn add_alternative_allele<I>(
        mut self,
        id: I,
        alternative_allele: Map<AlternativeAllele>,
    ) -> Self
    where
        I: Into<String>,
    {
        self.alternative_alleles
            .insert(id.into(), alternative_allele);

        self
    }

    /// Adds a contig record (`contig`).
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{self as vcf, header::record::value::{map::Contig, Map}};
    ///
    /// let contig = Map::<Contig>::new();
    ///
    /// let header = vcf::Header::builder()
    ///     .add_contig("sq0", contig.clone())
    ///     .build();
    ///
    /// let contigs = header.contigs();
    /// assert_eq!(contigs.len(), 1);
    /// assert_eq!(&contigs[0], &contig);
    /// ```
    pub fn add_contig<I>(mut self, id: I, contig: Map<Contig>) -> Self
    where
        I: Into<String>,
    {
        self.contigs.insert(id.into(), contig);
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
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::{value::Collection, Value},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .insert("fileDate".parse()?, Value::from("20200709"))?
    ///     .build();
    ///
    /// assert_eq!(
    ///     header.get("fileDate"),
    ///     Some(&Collection::Unstructured(vec![String::from("20200709")]))
    /// );
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn insert(
        mut self,
        key: record::key::Other,
        value: record::Value,
    ) -> Result<Self, super::record::value::collection::AddError> {
        let collection = self
            .other_records
            .entry(key)
            .or_insert_with(|| match value {
                record::Value::String(_) => record::value::Collection::Unstructured(Vec::new()),
                record::Value::Map(..) => record::value::Collection::Structured(IndexMap::new()),
            });

        collection.add(value)?;

        Ok(self)
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
            contigs: self.contigs,
            sample_names: self.sample_names,
            other_records: self.other_records,
            string_maps: StringMaps::default(),
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
        assert!(header.contigs().is_empty());
        assert!(header.sample_names().is_empty());
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header,
            variant::{
                record::info::field::key as info_key, record_buf::samples::keys::key as format_key,
            },
        };

        let (key, value) = (
            "fileDate".parse::<header::record::key::Other>()?,
            header::record::Value::from("20200709"),
        );

        let header = Builder::default()
            .set_file_format(FileFormat::new(4, 3))
            .add_info(
                info_key::SAMPLES_WITH_DATA_COUNT,
                Map::<Info>::from(info_key::SAMPLES_WITH_DATA_COUNT),
            )
            .add_filter("q10", Map::<Filter>::new("Quality below 10"))
            .add_format(
                format_key::GENOTYPE,
                Map::<Format>::from(format_key::GENOTYPE),
            )
            .add_alternative_allele("DEL", Map::<AlternativeAllele>::new("Deletion"))
            .add_contig("sq0", Map::<Contig>::new())
            .add_contig("sq1", Map::<Contig>::new())
            .add_sample_name("sample0")
            .insert(key.clone(), value.clone())?
            .insert(key.clone(), value)?
            .build();

        assert_eq!(header.file_format(), FileFormat::new(4, 3));
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.contigs().len(), 2);
        assert_eq!(header.get(&key).map(|collection| collection.len()), Some(2));

        Ok(())
    }
}
