use std::{
    ops::Deref,
    str::{FromStr, Lines},
};

use indexmap::IndexSet;
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
pub struct StringMap(IndexSet<String>);

impl StringMap {
    fn insert(&mut self, value: String) {
        self.0.insert(value);
    }

    fn insert_full(&mut self, value: String) -> (usize, bool) {
        self.0.insert_full(value)
    }
}

impl Default for StringMap {
    fn default() -> Self {
        // § 6.2.1 Dictionary of strings (2021-01-13): "Note that 'PASS' is always implicitly
        // encoded as the first entry in the header dictionary."
        Self([Filter::pass().id().into()].into_iter().collect())
    }
}

impl Deref for StringMap {
    type Target = IndexSet<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
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

            let (idx, j) = match record.key() {
                Key::Filter => {
                    let filter = Filter::try_from(record).map_err(ParseError::InvalidFilter)?;
                    let (j, _) = string_map.insert_full(filter.id().into());
                    (filter.idx(), j)
                }
                Key::Format => {
                    let format = Format::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidFormat)?;
                    let (j, _) = string_map.insert_full(format.id().as_ref().into());
                    (format.idx(), j)
                }
                Key::Info => {
                    let info = Info::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidInfo)?;
                    let (j, _) = string_map.insert_full(info.id().as_ref().into());
                    (info.idx(), j)
                }
                _ => continue,
            };

            if let Some(i) = idx {
                if i != j {
                    return Err(ParseError::StringMapPositionMismatch(i, j));
                }
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
        assert_eq!(
            StringMap::default(),
            StringMap([String::from("PASS")].into_iter().collect())
        );
    }

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
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

        let actual: StringMap = s.parse()?;
        let expected = StringMap(
            [
                String::from("PASS"),
                String::from("NS"),
                String::from("DP"),
                String::from("q10"),
                String::from("GT"),
            ]
            .into_iter()
            .collect(),
        );

        assert!(actual.iter().eq(expected.iter()));

        Ok(())
    }

    #[test]
    fn test_from_str_with_a_position_mismatch() {
        let s = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=8>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMap>(),
            Err(ParseError::StringMapPositionMismatch(8, 0))
        );

        let s = r#"##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=1>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMap>(),
            Err(ParseError::StringMapPositionMismatch(2, 1))
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
        let expected = StringMap(
            [
                String::from("PASS"),
                String::from("NS"),
                String::from("DP"),
                String::from("q10"),
                String::from("GT"),
            ]
            .into_iter()
            .collect(),
        );

        assert!(actual.iter().eq(expected.iter()));
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
