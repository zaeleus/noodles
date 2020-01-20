use std::{error, fmt};

use csv::StringRecord;

use super::{Attributes, Strand};

static EMPTY_VALUE: &str = ".";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Header {
    SeqName,
    Source,
    Feature,
    Start,
    End,
    Score,
    Strand,
    Frame,
    Attributes,
}

#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    Missing(Header),
    Empty(Header),
    Parse(Header, String),
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub type Result<T> = std::result::Result<T, Error>;

pub struct Record(StringRecord);

impl Record {
    pub fn new(inner: StringRecord) -> Self {
        Self(inner)
    }

    pub fn seq_name(&self) -> self::Result<&str> {
        self.parse(Header::SeqName)
    }

    pub fn source(&self) -> self::Result<&str> {
        self.parse(Header::Source)
    }

    pub fn feature(&self) -> self::Result<&str> {
        self.parse(Header::Feature)
    }

    pub fn start(&self) -> self::Result<u64> {
        self.parse_u64(Header::Start)
    }

    pub fn end(&self) -> self::Result<u64> {
        self.parse_u64(Header::End)
    }

    pub fn score(&self) -> self::Result<f64> {
        self.parse_f64(Header::Score)
    }

    pub fn strand(&self) -> self::Result<Strand> {
        self.parse(Header::Strand)
            .and_then(|s| s.parse().map_err(|e| Error::Parse(Header::Strand, e)))
    }

    pub fn frame(&self) -> self::Result<u8> {
        self.parse_u8(Header::Frame)
    }

    pub fn attributes(&self) -> self::Result<Attributes> {
        self.parse(Header::Attributes).map(Attributes::new)
    }

    pub fn into_inner(self) -> StringRecord {
        self.0
    }

    fn parse(&self, header: Header) -> self::Result<&str> {
        self.0
            .get(header as usize)
            .ok_or_else(|| Error::Missing(header))
            .and_then(|s| validate_empty(header, s))
    }

    fn parse_u64(&self, header: Header) -> self::Result<u64> {
        self.parse(header).and_then(|s| {
            s.parse()
                .map_err(|e| Error::Parse(header, format!("{}", e)))
        })
    }

    fn parse_f64(&self, header: Header) -> self::Result<f64> {
        self.parse(header).and_then(|s| {
            s.parse()
                .map_err(|e| Error::Parse(header, format!("{}", e)))
        })
    }

    fn parse_u8(&self, header: Header) -> self::Result<u8> {
        self.parse(header).and_then(|s| {
            s.parse()
                .map_err(|e| Error::Parse(header, format!("{}", e)))
        })
    }
}

fn validate_empty(header: Header, s: &str) -> self::Result<&str> {
    if is_empty(s) {
        Err(Error::Empty(header))
    } else {
        Ok(s)
    }
}

fn is_empty(s: &str) -> bool {
    s == EMPTY_VALUE
}

#[cfg(test)]
mod tests {
    use crate::Strand;

    use super::*;

    fn build_string_record() -> StringRecord {
        StringRecord::from(vec![
            "chr1",
            "HAVANA",
            "gene",
            "11869",
            "14409",
            "0.448",
            "+",
            "1",
            r#"gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"#,
        ])
    }

    fn build_empty_string_record() -> StringRecord {
        StringRecord::from(vec![".", ".", ".", ".", ".", ".", ".", "."])
    }

    #[test]
    fn test_seq_name() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.seq_name(), Ok("chr1"));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::SeqName);
        assert_eq!(record.seq_name(), Err(err));
    }

    #[test]
    fn test_source() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.source(), Ok("HAVANA"));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Source);
        assert_eq!(record.source(), Err(err));
    }

    #[test]
    fn test_feature() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.feature(), Ok("gene"));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Feature);
        assert_eq!(record.feature(), Err(err));
    }

    #[test]
    fn test_start() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.start(), Ok(11869));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Start);
        assert_eq!(record.start(), Err(err));
    }

    #[test]
    fn test_end() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.end(), Ok(14409));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::End);
        assert_eq!(record.end(), Err(err));
    }

    #[test]
    fn test_score() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.score(), Ok(0.448));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Score);
        assert_eq!(record.score(), Err(err));
    }

    #[test]
    fn test_strand() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.strand(), Ok(Strand::Forward));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Strand);
        assert_eq!(record.strand(), Err(err));
    }

    #[test]
    fn test_frame() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.frame(), Ok(1));

        let r = build_empty_string_record();
        let record = Record::new(r);
        let err = Error::Empty(Header::Frame);
        assert_eq!(record.frame(), Err(err));
    }

    #[test]
    fn test_attributes() {
        let r = build_string_record();
        let record = Record::new(r);
        assert!(record.attributes().is_ok());

        let r = build_empty_string_record();
        let record = Record::new(r);
        assert!(record.attributes().is_err());
    }

    #[test]
    fn test_validate_empty() {
        assert_eq!(
            validate_empty(Header::SeqName, "."),
            Err(Error::Empty(Header::SeqName))
        );
        assert_eq!(validate_empty(Header::Feature, "gene"), Ok("gene"));
        assert_eq!(validate_empty(Header::Start, "11869"), Ok("11869"));
        assert_eq!(validate_empty(Header::Score, "0.448"), Ok("0.448"));
        assert_eq!(validate_empty(Header::Strand, "+"), Ok("+"));
    }

    #[test]
    fn test_is_empty() {
        assert!(is_empty("."));
        assert!(!is_empty("gene"));
        assert!(!is_empty("11869"));
        assert!(!is_empty("0.448"));
        assert!(!is_empty("+"));
        assert!(!is_empty(r#"gene_id "ENSG00000223972.5";"#));
    }
}
