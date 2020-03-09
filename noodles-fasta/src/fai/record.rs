use csv::StringRecord;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Header {
    Name,
    Length,
    Offset,
    LineBases,
    LineWidth,
}

#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    Missing(Header),
    Parse(Header, String),
}

pub type Result<T> = std::result::Result<T, Error>;

pub struct Record(StringRecord);

impl Record {
    pub fn new(inner: StringRecord) -> Self {
        Self(inner)
    }

    pub fn name(&self) -> self::Result<&str> {
        self.parse(Header::Name)
    }

    pub fn length(&self) -> self::Result<u64> {
        self.parse_u64(Header::Length)
    }

    pub fn offset(&self) -> self::Result<u64> {
        self.parse_u64(Header::Offset)
    }

    pub fn line_bases(&self) -> self::Result<u64> {
        self.parse_u64(Header::LineBases)
    }

    pub fn line_width(&self) -> self::Result<u64> {
        self.parse_u64(Header::LineWidth)
    }

    fn parse(&self, header: Header) -> self::Result<&str> {
        self.0
            .get(header as usize)
            .ok_or_else(|| Error::Missing(header))
    }

    fn parse_u64(&self, header: Header) -> self::Result<u64> {
        self.parse(header).and_then(|s| {
            s.parse()
                .map_err(|e| Error::Parse(header, format!("{}", e)))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_string_record() -> StringRecord {
        StringRecord::from(vec!["chr1", "248956422", "112", "70", "71"])
    }

    #[test]
    fn test_name() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.name(), Ok("chr1"));
    }

    #[test]
    fn test_length() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.length(), Ok(248956422));
    }

    #[test]
    fn test_offset() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.offset(), Ok(112));
    }

    #[test]
    fn test_line_bases() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.line_bases(), Ok(70));
    }

    #[test]
    fn test_line_width() {
        let r = build_string_record();
        let record = Record::new(r);
        assert_eq!(record.line_width(), Ok(71));
    }
}
