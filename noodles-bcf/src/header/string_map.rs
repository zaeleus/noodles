use std::{
    collections::HashMap,
    mem,
    str::{FromStr, Lines},
};

use noodles_vcf::{
    self as vcf,
    header::{Filter, Format, Info, ParseError, Record},
};

/// An indexed map of VCF strings.
///
/// This is also called a dictionary of strings.
///
/// See § 6.2.1 Dictionary of strings (2021-05-13).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StringMap {
    indices: HashMap<String, usize>,
    entries: Vec<Option<String>>,
}

impl StringMap {
    /// Returns an entry by index.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::header::StringMap;
    /// let string_map = StringMap::default();
    /// assert_eq!(string_map.get_index(0), Some("PASS"));
    /// assert!(string_map.get_index(1).is_none());
    /// ```
    pub fn get_index(&self, i: usize) -> Option<&str> {
        self.entries.get(i).and_then(|entry| entry.as_deref())
    }

    /// Returns the index of the entry of the given value.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::header::StringMap;
    /// let string_map = StringMap::default();
    /// assert_eq!(string_map.get_index_of("PASS"), Some(0));
    /// assert!(string_map.get_index_of("DP").is_none());
    /// ```
    pub fn get_index_of(&self, value: &str) -> Option<usize> {
        self.indices.get(value).copied()
    }

    fn get_full(&self, value: &str) -> Option<(usize, &str)> {
        self.get_index_of(value)
            .and_then(|i| self.get_index(i).map(|entry| (i, entry)))
    }

    fn insert(&mut self, value: String) -> Option<String> {
        self.insert_full(value).1
    }

    fn insert_full(&mut self, value: String) -> (usize, Option<String>) {
        match self.get_index_of(&value) {
            Some(i) => {
                let entry = mem::replace(&mut self.entries[i], Some(value));
                (i, entry)
            }
            None => {
                let i = self.push(value);
                (i, None)
            }
        }
    }

    fn insert_at(&mut self, i: usize, value: String) -> Option<String> {
        if i >= self.entries.len() {
            self.entries.resize(i + 1, None);
        }

        self.indices.insert(value.clone(), i);
        mem::replace(&mut self.entries[i], Some(value))
    }

    fn push(&mut self, value: String) -> usize {
        let i = self.entries.len();

        self.indices.insert(value.clone(), i);
        self.entries.push(Some(value));

        i
    }
}

impl Default for StringMap {
    fn default() -> Self {
        // § 6.2.1 Dictionary of strings (2021-01-13): "Note that 'PASS' is always implicitly
        // encoded as the first entry in the header dictionary."
        let pass = Filter::pass().id().to_string();

        Self {
            indices: [(pass.clone(), 0)].into_iter().collect(),
            entries: vec![Some(pass)],
        }
    }
}

impl FromStr for StringMap {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use vcf::header::record::Key;

        let mut string_map = StringMap::default();

        let mut lines = s.lines();
        let file_format = parse_file_format(&mut lines)?;

        for line in &mut lines {
            if line.starts_with("#CHROM") {
                break;
            }

            let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;

            let (id, idx) = match record.key() {
                Key::Filter => {
                    let filter = Filter::try_from(record).map_err(ParseError::InvalidFilter)?;
                    (filter.id().to_string(), filter.idx())
                }
                Key::Format => {
                    let format = Format::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidFormat)?;
                    (format.id().as_ref().to_string(), format.idx())
                }
                Key::Info => {
                    let info = Info::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidInfo)?;
                    (info.id().as_ref().to_string(), info.idx())
                }
                _ => continue,
            };

            if let Some(i) = idx {
                if let Some((j, entry)) = string_map.get_full(&id) {
                    let actual = (i, id);
                    let expected = (j, entry.into());

                    if actual != expected {
                        return Err(ParseError::StringMapPositionMismatch(actual, expected));
                    }
                } else {
                    string_map.insert_at(i, id);
                }
            } else {
                string_map.insert(id);
            }
        }

        Ok(string_map)
    }
}

fn parse_file_format(lines: &mut Lines<'_>) -> Result<vcf::header::FileFormat, ParseError> {
    use vcf::header::record::{Key, Value};

    let record: Record = lines
        .next()
        .ok_or(ParseError::MissingFileFormat)
        .and_then(|line| line.parse().map_err(ParseError::InvalidRecord))?;

    if record.key() == &Key::FileFormat {
        match record.value() {
            Value::String(value) => value.parse().map_err(ParseError::InvalidFileFormat),
            _ => Err(ParseError::InvalidRecordValue),
        }
    } else {
        Err(ParseError::MissingFileFormat)
    }
}

impl From<&vcf::Header> for StringMap {
    fn from(header: &vcf::Header) -> Self {
        let mut string_map = StringMap::default();

        for info in header.infos().values() {
            string_map.insert(info.id().as_ref().into());
        }

        for filter in header.filters().values() {
            string_map.insert(filter.id().into());
        }

        for format in header.formats().values() {
            string_map.insert(format.id().as_ref().into());
        }

        string_map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let string_map = StringMap::default();

        assert_eq!(string_map.indices.len(), 1);
        assert_eq!(string_map.indices.get("PASS"), Some(&0));

        assert_eq!(string_map.entries.len(), 1);
        assert_eq!(string_map.entries[0], Some(String::from("PASS")));
    }

    #[test]
    fn test_from_str() {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20210412
##contig=<ID=sq0,length=8>
##contig=<ID=sq1,length=13>
##contig=<ID=sq2,length=21>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data",IDX=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=2>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Quality below 10",IDX=3>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=4>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=2>
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let indices = [
            (String::from("PASS"), 0),
            (String::from("NS"), 1),
            (String::from("DP"), 2),
            (String::from("q10"), 3),
            (String::from("GT"), 4),
        ]
        .into_iter()
        .collect();
        let entries = vec![
            Some(String::from("PASS")),
            Some(String::from("NS")),
            Some(String::from("DP")),
            Some(String::from("q10")),
            Some(String::from("GT")),
        ];
        let expected = StringMap { indices, entries };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_mixed_positions() {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20210412
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data",IDX=1>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Quality below 10",IDX=3>
##FILTER=<ID=q15,Description="Quality below 15",IDX=4>
##FILTER=<ID=q20,Description="Quality below 20">
##FILTER=<ID=NS,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let indices = [
            (String::from("PASS"), 0),
            (String::from("NS"), 1),
            (String::from("q10"), 3),
            (String::from("q15"), 4),
            (String::from("q20"), 5),
        ]
        .into_iter()
        .collect();
        let entries = vec![
            Some(String::from("PASS")),
            Some(String::from("NS")),
            None,
            Some(String::from("q10")),
            Some(String::from("q15")),
            Some(String::from("q20")),
        ];
        let expected = StringMap { indices, entries };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_a_position_mismatch() {
        let s = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=8>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMap>(),
            Err(ParseError::StringMapPositionMismatch(
                (8, String::from("PASS")),
                (0, String::from("PASS"))
            ))
        );

        let s = r#"##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=1>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMap>(),
            Err(ParseError::StringMapPositionMismatch(
                (2, String::from("DP")),
                (1, String::from("DP"))
            ))
        );
    }

    #[test]
    fn test_vcf_header_for_string_map() {
        use vcf::{
            header::{AlternativeAllele, Contig},
            record::{
                alternate_bases::allele::{
                    symbol::{structural_variant::Type, StructuralVariant},
                    Symbol,
                },
                genotypes::genotype::field::Key as GenotypeKey,
                info::field::Key as InfoKey,
            },
        };

        let header = vcf::Header::builder()
            .add_contig(Contig::new("sq0"))
            .add_contig(Contig::new("sq1"))
            .add_contig(Contig::new("sq2"))
            .add_info(Info::from(InfoKey::SamplesWithDataCount))
            .add_info(Info::from(vcf::record::info::field::Key::TotalDepth))
            .add_filter(Filter::pass())
            .add_filter(Filter::new("q10", "Quality below 10"))
            .add_format(Format::from(GenotypeKey::Genotype))
            .add_format(Format::from(GenotypeKey::ReadDepth))
            .add_alternative_allele(AlternativeAllele::new(
                Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
                "Deletion",
            ))
            .build();

        let actual = StringMap::from(&header);

        let indices = [
            (String::from("PASS"), 0),
            (String::from("NS"), 1),
            (String::from("DP"), 2),
            (String::from("q10"), 3),
            (String::from("GT"), 4),
        ]
        .into_iter()
        .collect();
        let entries = vec![
            Some(String::from("PASS")),
            Some(String::from("NS")),
            Some(String::from("DP")),
            Some(String::from("q10")),
            Some(String::from("GT")),
        ];
        let expected = StringMap { indices, entries };

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_parse_file_format() {
        use vcf::header::FileFormat;

        let s = "##fileformat=VCFv4.3\n";
        let mut lines = s.lines();
        assert_eq!(parse_file_format(&mut lines), Ok(FileFormat::new(4, 3)));

        let s = "";
        let mut lines = s.lines();
        assert_eq!(
            parse_file_format(&mut lines),
            Err(ParseError::MissingFileFormat)
        );

        let s = "fileformat=VCFv4.3";
        let mut lines = s.lines();
        assert!(matches!(
            parse_file_format(&mut lines),
            Err(ParseError::InvalidRecord(_))
        ));

        let s = "##fileformat=VCF43\n";
        let mut lines = s.lines();
        assert!(matches!(
            parse_file_format(&mut lines),
            Err(ParseError::InvalidFileFormat(_))
        ));

        let s = "##fileDate=20211119";
        let mut lines = s.lines();
        assert_eq!(
            parse_file_format(&mut lines),
            Err(ParseError::MissingFileFormat)
        );
    }
}
