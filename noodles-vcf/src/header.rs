mod filter;
mod info;
mod record;

pub use self::{filter::Filter, info::Info};

use std::{convert::TryFrom, str::FromStr};

use self::record::Record;

#[derive(Debug, Default)]
pub struct Header {
    infos: Vec<Info>,
    filters: Vec<Filter>,
}

impl Header {
    pub fn infos(&self) -> &[Info] {
        &self.infos
    }

    pub fn filters(&self) -> &[Filter] {
        &self.filters
    }
}

#[derive(Debug)]
pub enum ParseError {
    InvalidRecord(record::ParseError),
    InvalidInfo(info::ParseError),
    InvalidFilter(filter::ParseError),
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut header = Header::default();

        for line in s.lines() {
            println!("{}", line);
            let record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record {
                Record::Info(fields) => {
                    let info = Info::try_from(&fields[..]).map_err(ParseError::InvalidInfo)?;
                    header.infos.push(info);
                }
                Record::Filter(fields) => {
                    let filter =
                        Filter::try_from(&fields[..]).map_err(ParseError::InvalidFilter)?;
                    header.filters.push(filter);
                }
            }
        }

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##FILTER=<ID=q10,Description="Quality below 10">
"#;

        let header: Header = s.parse()?;

        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);

        Ok(())
    }
}
