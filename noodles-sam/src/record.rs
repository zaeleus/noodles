mod builder;
mod field;
mod quality_scores;
mod read_name;
mod reference_sequence_name;
mod sequence;

pub use self::{
    builder::Builder, field::Field, quality_scores::QualityScores, read_name::ReadName,
    reference_sequence_name::ReferenceSequenceName, sequence::Sequence,
};

use std::{error, fmt, str::FromStr};

use super::{Cigar, Data, Flags, MappingQuality};

pub(crate) const NULL_FIELD: &str = "*";
const FIELD_DELIMITER: char = '\t';
const MAX_FIELDS: usize = 12;

#[derive(Debug)]
pub struct Record {
    qname: ReadName,
    flag: Flags,
    rname: ReferenceSequenceName,
    pos: u32,
    mapq: MappingQuality,
    cigar: Cigar,
    rnext: String,
    pnext: u32,
    tlen: i32,
    seq: Sequence,
    qual: QualityScores,
    data: Data,
}

impl Record {
    pub fn builder() -> Builder {
        Builder::new()
    }

    pub fn name(&self) -> &ReadName {
        &self.qname
    }

    pub fn flags(&self) -> Flags {
        self.flag
    }

    pub fn reference_sequence_name(&self) -> &ReferenceSequenceName {
        &self.rname
    }

    pub fn position(&self) -> u32 {
        self.pos
    }

    pub fn mapping_quality(&self) -> MappingQuality {
        self.mapq
    }

    pub fn cigar(&self) -> &Cigar {
        &self.cigar
    }

    pub fn mate_reference_sequence_name(&self) -> &str {
        &self.rnext
    }

    pub fn mate_position(&self) -> u32 {
        self.pnext
    }

    pub fn template_len(&self) -> i32 {
        self.tlen
    }

    pub fn sequence(&self) -> &Sequence {
        &self.seq
    }

    pub fn quality_scores(&self) -> &QualityScores {
        &self.qual
    }

    pub fn data(&self) -> &Data {
        &self.data
    }
}

impl Default for Record {
    fn default() -> Self {
        Builder::new().build()
    }
}

#[derive(Debug)]
pub enum ParseError {
    Missing(Field),
    Invalid(Field, String),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(field) => write!(f, "missing field: {}", field.name()),
            Self::Invalid(field, message) => {
                write!(f, "invalid {} field: {}", field.name(), message)
            }
        }
    }
}

impl FromStr for Record {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.splitn(MAX_FIELDS, FIELD_DELIMITER);

        let qname = parse_string(&mut fields, Field::Name).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::Name, format!("{}", e)))
        })?;

        let flag = parse_u16(&mut fields, Field::Flags).map(Flags::from)?;

        let rname = parse_string(&mut fields, Field::ReferenceSequenceName).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::ReferenceSequenceName, format!("{}", e)))
        })?;

        let pos = parse_u32(&mut fields, Field::Position)?;
        let mapq = parse_u8(&mut fields, Field::MappingQuality).map(MappingQuality::from)?;

        let cigar = parse_string(&mut fields, Field::Cigar).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::Cigar, format!("{}", e)))
        })?;

        let rnext = parse_string(&mut fields, Field::MateReferenceSequenceName)?;
        let pnext = parse_u32(&mut fields, Field::MatePosition)?;
        let tlen = parse_i32(&mut fields, Field::TemplateLength)?;

        let seq = parse_string(&mut fields, Field::Sequence).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::Sequence, format!("{}", e)))
        })?;

        let qual = parse_string(&mut fields, Field::QualityScores).and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(Field::QualityScores, format!("{}", e)))
        })?;

        let data = match fields.next() {
            Some(s) => s
                .parse()
                .map_err(|e| ParseError::Invalid(Field::Data, format!("{}", e)))?,
            None => Data::default(),
        };

        Ok(Record {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            data,
        })
    }
}

fn parse_string<'a, I>(fields: &mut I, field: Field) -> Result<String, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .map(|s| s.into())
}

fn parse_u8<'a, I>(fields: &mut I, field: Field) -> Result<u8, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(field, format!("{}", e)))
        })
}

fn parse_u16<'a, I>(fields: &mut I, field: Field) -> Result<u16, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(field, format!("{}", e)))
        })
}

fn parse_i32<'a, I>(fields: &mut I, field: Field) -> Result<i32, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(field, format!("{}", e)))
        })
}

fn parse_u32<'a, I>(fields: &mut I, field: Field) -> Result<u32, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or_else(|| ParseError::Missing(field))
        .and_then(|s| {
            s.parse()
                .map_err(|e| ParseError::Invalid(field, format!("{}", e)))
        })
}
