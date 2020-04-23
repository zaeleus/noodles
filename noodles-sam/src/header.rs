#[allow(clippy::module_inception)]
pub mod header;
mod program;
mod read_group;
mod record;
mod reference_sequence;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use indexmap::IndexMap;

pub use self::{program::Program, read_group::ReadGroup, reference_sequence::ReferenceSequence};

pub use self::record::Record;

pub type ReferenceSequences = IndexMap<String, ReferenceSequence>;

#[derive(Debug, Default)]
pub struct Header {
    header: header::Header,
    reference_sequences: ReferenceSequences,
    read_groups: Vec<ReadGroup>,
    programs: Vec<Program>,
    comments: Vec<String>,
}

impl Header {
    pub fn header(&self) -> &header::Header {
        &self.header
    }

    pub fn reference_sequences(&self) -> &ReferenceSequences {
        &self.reference_sequences
    }

    pub fn read_groups(&self) -> &[ReadGroup] {
        &self.read_groups
    }

    pub fn programs(&self) -> &[Program] {
        &self.programs
    }

    pub fn comments(&self) -> &[String] {
        &self.comments
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{}", self.header)?;

        for reference_sequence in self.reference_sequences.values() {
            writeln!(f, "{}", reference_sequence)?;
        }

        for read_group in &self.read_groups {
            writeln!(f, "{}", read_group)?;
        }

        for program in &self.programs {
            writeln!(f, "{}", program)?;
        }

        for comment in &self.comments {
            writeln!(f, "{}\t{}", record::Kind::Comment, comment)?;
        }

        Ok(())
    }
}

#[derive(Debug)]
pub enum ParseError {
    InvalidRecord(record::ParseError),
    InvalidHeader(header::ParseError),
    InvalidReferenceSequence(reference_sequence::ParseError),
    InvalidReadGroup(read_group::ParseError),
    InvalidProgram(program::ParseError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord(e) => write!(f, "{}", e),
            Self::InvalidHeader(e) => write!(f, "{}", e),
            Self::InvalidReferenceSequence(e) => write!(f, "{}", e),
            Self::InvalidReadGroup(e) => write!(f, "{}", e),
            Self::InvalidProgram(e) => write!(f, "{}", e),
        }
    }
}

impl FromStr for Header {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut header = Header::default();

        for line in s.lines() {
            let record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record {
                Record::Header(fields) => {
                    header.header =
                        header::Header::try_from(&fields[..]).map_err(ParseError::InvalidHeader)?;
                }
                Record::ReferenceSequence(fields) => {
                    let reference_sequence = ReferenceSequence::try_from(&fields[..])
                        .map_err(ParseError::InvalidReferenceSequence)?;

                    let name = reference_sequence.name();

                    header
                        .reference_sequences
                        .insert(name.into(), reference_sequence);
                }
                Record::ReadGroup(fields) => {
                    let read_group =
                        ReadGroup::try_from(&fields[..]).map_err(ParseError::InvalidReadGroup)?;
                    header.read_groups.push(read_group);
                }
                Record::Program(fields) => {
                    let program =
                        Program::try_from(&fields[..]).map_err(ParseError::InvalidProgram)?;
                    header.programs.push(program);
                }
                Record::Comment(comment) => {
                    header.comments.push(comment);
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
    fn test_fmt() {
        let reference_sequences = vec![
            (
                String::from("sq0"),
                ReferenceSequence::new(String::from("sq0"), 8),
            ),
            (
                String::from("sq1"),
                ReferenceSequence::new(String::from("sq1"), 13),
            ),
        ]
        .into_iter()
        .collect();

        let header = Header {
            header: header::Header::new(String::from("1.6")),
            reference_sequences,
            read_groups: vec![
                ReadGroup::new(String::from("rg0")),
                ReadGroup::new(String::from("rg1")),
            ],
            programs: vec![
                Program::new(String::from("pg0")),
                Program::new(String::from("pg1")),
            ],
            comments: vec![String::from("noodles"), String::from("sam")],
        };

        let actual = format!("{}", header);
        let expected = "\
@HD\tVN:1.6
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@RG\tID:rg1
@PG\tID:pg0
@PG\tID:pg1
@CO\tnoodles
@CO\tsam
";

        assert_eq!(actual, expected);
    }

    #[test]
    fn test_from_str() {
        let raw_header = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:1
@SQ\tSN:sq1\tLN:2
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tnoodles_sam::header::tests::test_from_str
";

        let header: Header = raw_header.parse().unwrap();

        assert_eq!(header.reference_sequences().len(), 2);

        assert_eq!(header.read_groups().len(), 1);

        assert_eq!(header.programs().len(), 1);

        assert_eq!(header.comments.len(), 1);
        assert_eq!(
            &header.comments[0],
            "noodles_sam::header::tests::test_from_str"
        );
    }
}
