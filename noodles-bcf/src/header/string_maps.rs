//! An indexed map of VCF strings.

mod string_map;

use std::str::{FromStr, Lines};

use noodles_vcf::{
    self as vcf,
    header::{Contig, Filter, Format, Info, ParseError, Record},
};

pub use self::string_map::StringMap;

/// An indexed map of VCF strings (FILTER, FORMAT, and INFO).
pub type StringStringMap = StringMap;

/// An indexed map of VCF contig names.
pub type ContigStringMap = StringMap;

/// An indexed map of VCF strings.
///
/// This includes both the dictionary of strings and dictionary of contigs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StringMaps {
    string_string_map: StringStringMap,
    contig_string_map: ContigStringMap,
}

impl StringMaps {
    /// Returns an indexed map of VCF strings (FILTER, FORMAT, and INFO).
    ///
    /// The filter ID "PASS" is always the first entry in the string string map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::{self as bcf, header::StringMaps};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Key as FormatKey, info::Key as InfoKey, Contig, Filter, Format, Info},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(Info::from(InfoKey::TotalDepth))
    ///     .add_filter(Filter::new("q10", "Quality below 10"))
    ///     .add_format(Format::from(FormatKey::ReadDepth))
    ///     .add_contig(Contig::new("sq0"))
    ///     .build();
    ///
    /// let string_maps = StringMaps::from(&header);
    /// let string_string_map = string_maps.strings();
    ///
    /// assert_eq!(string_string_map.get_index(0), Some("PASS"));
    /// assert_eq!(string_string_map.get_index(1), Some("DP"));
    /// assert_eq!(string_string_map.get_index(2), Some("q10"));
    /// assert!(string_string_map.get_index(3).is_none());
    /// # Ok::<_, vcf::header::ParseError>(())
    /// ```
    pub fn strings(&self) -> &StringStringMap {
        &self.string_string_map
    }

    fn strings_mut(&mut self) -> &mut StringStringMap {
        &mut self.string_string_map
    }

    /// Returns an indexed map of contig names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::{self as bcf, header::StringMaps};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{format::Key as FormatKey, info::Key as InfoKey, Contig, Filter, Format, Info},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(Info::from(InfoKey::TotalDepth))
    ///     .add_filter(Filter::new("q10", "Quality below 10"))
    ///     .add_format(Format::from(FormatKey::ReadDepth))
    ///     .add_contig(Contig::new("sq0"))
    ///     .build();
    ///
    /// let string_maps = StringMaps::from(&header);
    /// let contig_string_map = string_maps.contigs();
    ///
    /// assert_eq!(contig_string_map.get_index(0), Some("sq0"));
    /// assert!(contig_string_map.get_index(1).is_none());
    /// # Ok::<_, vcf::header::ParseError>(())
    /// ```
    pub fn contigs(&self) -> &ContigStringMap {
        &self.contig_string_map
    }

    fn contigs_mut(&mut self) -> &mut ContigStringMap {
        &mut self.contig_string_map
    }
}

impl Default for StringMaps {
    fn default() -> Self {
        // § 6.2.1 Dictionary of strings (2021-01-13): "Note that 'PASS' is always implicitly
        // encoded as the first entry in the header dictionary."
        let mut string_string_map = StringMap::default();
        let pass = Filter::pass().id().to_string();
        string_string_map.insert(pass);

        let contig_string_map = StringMap::default();

        Self {
            string_string_map,
            contig_string_map,
        }
    }
}

impl FromStr for StringMaps {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use vcf::header::record::Key;

        let mut string_maps = Self::default();

        let mut lines = s.lines();
        let file_format = parse_file_format(&mut lines)?;

        for line in &mut lines {
            if line.starts_with("#CHROM") {
                break;
            }

            let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record.key() {
                Key::Contig => {
                    let contig = Contig::try_from(record).map_err(ParseError::InvalidContig)?;
                    insert(string_maps.contigs_mut(), contig.id(), contig.idx())?;
                }
                Key::Filter => {
                    let filter = Filter::try_from(record).map_err(ParseError::InvalidFilter)?;
                    insert(string_maps.strings_mut(), filter.id(), filter.idx())?;
                }
                Key::Format => {
                    let format = Format::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidFormat)?;

                    insert(
                        string_maps.strings_mut(),
                        format.id().as_ref(),
                        format.idx(),
                    )?;
                }
                Key::Info => {
                    let info = Info::try_from_record_file_format(record, file_format)
                        .map_err(ParseError::InvalidInfo)?;
                    insert(string_maps.strings_mut(), info.id().as_ref(), info.idx())?;
                }
                _ => {}
            }
        }

        Ok(string_maps)
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

fn insert(string_map: &mut StringMap, id: &str, idx: Option<usize>) -> Result<(), ParseError> {
    if let Some(i) = idx {
        if let Some((j, entry)) = string_map.get_full(id) {
            let actual = (i, id.into());
            let expected = (j, entry.into());

            if actual != expected {
                return Err(ParseError::StringMapPositionMismatch(actual, expected));
            }
        } else {
            string_map.insert_at(i, id.into());
        }
    } else {
        string_map.insert(id.into());
    }

    Ok(())
}

impl From<&vcf::Header> for StringMaps {
    fn from(header: &vcf::Header) -> Self {
        let mut string_maps = StringMaps::default();

        for contig in header.contigs().values() {
            string_maps.contigs_mut().insert(contig.id().into());
        }

        for info in header.infos().values() {
            string_maps.strings_mut().insert(info.id().as_ref().into());
        }

        for filter in header.filters().values() {
            string_maps.strings_mut().insert(filter.id().into());
        }

        for format in header.formats().values() {
            string_maps
                .strings_mut()
                .insert(format.id().as_ref().into());
        }

        string_maps
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let actual = StringMaps::default();

        let mut string_string_map = StringMap::default();
        string_string_map.insert("PASS".into());

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(actual, expected);
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

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("DP"), 2),
                (String::from("q10"), 3),
                (String::from("GT"), 4),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                Some(String::from("DP")),
                Some(String::from("q10")),
                Some(String::from("GT")),
            ],
        };

        let contig_string_map = StringMap {
            indices: [
                (String::from("sq0"), 0),
                (String::from("sq1"), 1),
                (String::from("sq2"), 2),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("sq0")),
                Some(String::from("sq1")),
                Some(String::from("sq2")),
            ],
        };

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

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

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("q10"), 3),
                (String::from("q15"), 4),
                (String::from("q20"), 5),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                None,
                Some(String::from("q10")),
                Some(String::from("q15")),
                Some(String::from("q20")),
            ],
        };

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_a_position_mismatch() {
        let s = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=8>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMaps>(),
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
            s.parse::<StringMaps>(),
            Err(ParseError::StringMapPositionMismatch(
                (2, String::from("DP")),
                (1, String::from("DP"))
            ))
        );
    }

    #[test]
    fn test_vcf_header_for_string_map() {
        use vcf::{
            header::{format::Key as FormatKey, info::Key as InfoKey, AlternativeAllele, Contig},
            record::alternate_bases::allele::{
                symbol::{structural_variant::Type, StructuralVariant},
                Symbol,
            },
        };

        let header = vcf::Header::builder()
            .add_contig(Contig::new("sq0"))
            .add_contig(Contig::new("sq1"))
            .add_contig(Contig::new("sq2"))
            .add_info(Info::from(InfoKey::SamplesWithDataCount))
            .add_info(Info::from(InfoKey::TotalDepth))
            .add_filter(Filter::pass())
            .add_filter(Filter::new("q10", "Quality below 10"))
            .add_format(Format::from(FormatKey::Genotype))
            .add_format(Format::from(FormatKey::ReadDepth))
            .add_alternative_allele(AlternativeAllele::new(
                Symbol::StructuralVariant(StructuralVariant::from(Type::Deletion)),
                "Deletion",
            ))
            .build();

        let actual = StringMaps::from(&header);

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("DP"), 2),
                (String::from("q10"), 3),
                (String::from("GT"), 4),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                Some(String::from("DP")),
                Some(String::from("q10")),
                Some(String::from("GT")),
            ],
        };

        let contig_string_map = StringMap {
            indices: [
                (String::from("sq0"), 0),
                (String::from("sq1"), 1),
                (String::from("sq2"), 2),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("sq0")),
                Some(String::from("sq1")),
                Some(String::from("sq2")),
            ],
        };

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

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
