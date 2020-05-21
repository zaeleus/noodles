mod builder;
#[allow(clippy::module_inception)]
pub mod header;
pub mod program;
pub mod read_group;
mod record;
pub mod reference_sequence;

use std::{convert::TryFrom, error, fmt, str::FromStr};

use indexmap::IndexMap;

pub use self::{
    builder::Builder, program::Program, read_group::ReadGroup,
    reference_sequence::ReferenceSequence,
};

pub use self::record::Record;

pub type ReferenceSequences = IndexMap<String, ReferenceSequence>;

#[derive(Debug, Default)]
pub struct Header {
    header: Option<header::Header>,
    reference_sequences: ReferenceSequences,
    read_groups: Vec<ReadGroup>,
    programs: Vec<Program>,
    comments: Vec<String>,
}

impl Header {
    pub fn builder() -> Builder {
        Builder::new()
    }

    pub fn header(&self) -> Option<&header::Header> {
        self.header.as_ref()
    }

    pub fn header_mut(&mut self) -> Option<&mut header::Header> {
        self.header.as_mut()
    }

    pub fn reference_sequences(&self) -> &ReferenceSequences {
        &self.reference_sequences
    }

    pub fn reference_sequences_mut(&mut self) -> &mut ReferenceSequences {
        &mut self.reference_sequences
    }

    pub fn read_groups(&self) -> &[ReadGroup] {
        &self.read_groups
    }

    pub fn read_groups_mut(&mut self) -> &mut Vec<ReadGroup> {
        &mut self.read_groups
    }

    pub fn programs(&self) -> &[Program] {
        &self.programs
    }

    pub fn programs_mut(&mut self) -> &mut Vec<Program> {
        &mut self.programs
    }

    pub fn comments(&self) -> &[String] {
        &self.comments
    }

    pub fn comments_mut(&mut self) -> &mut Vec<String> {
        &mut self.comments
    }
}

impl fmt::Display for Header {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(header) = self.header() {
            writeln!(f, "{}", header)?;
        }

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
    UnexpectedHeader,
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
            Self::UnexpectedHeader => f.write_str("unexpected @HD"),
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

        for (i, line) in s.lines().enumerate() {
            let record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record {
                Record::Header(fields) => {
                    if i == 0 {
                        header.header = Some(
                            header::Header::try_from(&fields[..])
                                .map_err(ParseError::InvalidHeader)?,
                        );
                    } else {
                        return Err(ParseError::UnexpectedHeader);
                    }
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
        let header = Header::builder()
            .set_header(header::Header::new(String::from("1.6")))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq0"), 8))
            .add_reference_sequence(ReferenceSequence::new(String::from("sq1"), 13))
            .add_read_group(ReadGroup::new(String::from("rg0")))
            .add_read_group(ReadGroup::new(String::from("rg1")))
            .add_program(Program::new(String::from("pg0")))
            .add_program(Program::new(String::from("pg1")))
            .add_comment("noodles")
            .add_comment("sam")
            .build();

        let actual = header.to_string();
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
    fn test_from_str() -> Result<(), ParseError> {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:sq0\tLN:1
@SQ\tSN:sq1\tLN:2
@RG\tID:rg0
@PG\tID:pg0\tPN:noodles
@CO\tnoodles_sam::header::tests::test_from_str
";

        let header: Header = s.parse()?;

        assert_eq!(header.reference_sequences().len(), 2);
        assert_eq!(header.read_groups().len(), 1);
        assert_eq!(header.programs().len(), 1);
        assert_eq!(header.comments.len(), 1);
        assert_eq!(
            &header.comments[0],
            "noodles_sam::header::tests::test_from_str"
        );

        Ok(())
    }

    #[test]
    fn test_from_str_with_empty_input() -> Result<(), ParseError> {
        let header: Header = "".parse()?;

        assert!(header.header().is_none());
        assert!(header.reference_sequences().is_empty());
        assert!(header.read_groups().is_empty());
        assert!(header.programs().is_empty());
        assert!(header.comments().is_empty());

        Ok(())
    }

    #[test]
    fn test_from_str_without_hd() -> Result<(), ParseError> {
        let header: Header = "@SQ\tSN:sq0\tLN:8\n".parse()?;
        assert!(header.header().is_none());
        assert_eq!(header.reference_sequences().len(), 1);
        Ok(())
    }

    #[test]
    fn test_from_str_with_multiple_hd() {
        let s = "\
@HD\tVN:1.6\tSO:coordinate
@HD\tVN:1.6\tSO:coordinate
";

        assert!(s.parse::<Header>().is_err());
    }
}
