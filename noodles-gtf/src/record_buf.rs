//! GTF record and fields.

pub mod attributes;
mod builder;
mod convert;
pub mod frame;
pub mod strand;

pub use self::{attributes::Attributes, builder::Builder, frame::Frame, strand::Strand};

use std::{error, fmt, num, str::FromStr};

use noodles_core::Position;

pub(crate) const MISSING_FIELD: &str = ".";

/// A GTF record buffer.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBuf {
    reference_sequence_name: String,
    source: String,
    ty: String,
    start: Position,
    end: Position,
    score: Option<f32>,
    strand: Option<Strand>,
    frame: Option<Frame>,
    attributes: Attributes,
}

impl RecordBuf {
    /// Returns a record builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let builder = gtf::RecordBuf::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the reference sequence name.
    ///
    /// This is also called the "seqname".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.reference_sequence_name(), ".");
    /// ```
    pub fn reference_sequence_name(&self) -> &str {
        &self.reference_sequence_name
    }

    /// Returns the source.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.source(), ".");
    /// ```
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Returns the feature type.
    ///
    /// This is also simply called "feature".
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.ty(), ".");
    /// ```
    pub fn ty(&self) -> &str {
        &self.ty
    }

    /// Returns the start position.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.start(), Position::MIN);
    /// ```
    pub fn start(&self) -> Position {
        self.start
    }

    /// Returns the end position.
    ///
    /// This position is 1-based, inclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert_eq!(record.end(), Position::MIN);
    /// ```
    pub fn end(&self) -> Position {
        self.end
    }

    /// Returns the confidence score.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.score().is_none());
    /// ```
    pub fn score(&self) -> Option<f32> {
        self.score
    }

    /// Returns the strand.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.strand().is_none());
    /// ```
    pub fn strand(&self) -> Option<Strand> {
        self.strand
    }

    /// Returns the frame.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.frame().is_none());
    /// ```
    pub fn frame(&self) -> Option<Frame> {
        self.frame
    }

    /// Returns the attributes.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_gtf as gtf;
    /// let record = gtf::RecordBuf::default();
    /// assert!(record.attributes().is_empty());
    /// ```
    pub fn attributes(&self) -> &Attributes {
        &self.attributes
    }
}

impl Default for RecordBuf {
    fn default() -> Self {
        Self::builder().build()
    }
}

/// An error returned when a raw GTF record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The source is missing.
    MissingSource,
    /// The type is missing.
    MissingType,
    /// The start is missing.
    MissingStart,
    /// The start is invalid.
    InvalidStart(num::ParseIntError),
    /// The end is missing.
    MissingEnd,
    /// The end is invalid.
    InvalidEnd(num::ParseIntError),
    /// The score is missing.
    MissingScore,
    /// The score is invalid.
    InvalidScore(num::ParseFloatError),
    /// The strand is missing.
    MissingStrand,
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
    /// The frame is missing.
    MissingFrame,
    /// The frame is invalid.
    InvalidFrame(frame::ParseError),
    /// The attributes are missing.
    MissingAttributes,
    /// The attributes are invalid.
    InvalidAttributes(attributes::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStart(e) | Self::InvalidEnd(e) => Some(e),
            Self::InvalidScore(e) => Some(e),
            Self::InvalidStrand(e) => Some(e),
            Self::InvalidFrame(e) => Some(e),
            Self::InvalidAttributes(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty input"),
            Self::MissingReferenceSequenceName => write!(f, "missing reference sequence name"),
            Self::MissingSource => write!(f, "missing source"),
            Self::MissingType => write!(f, "missing type"),
            Self::MissingStart => write!(f, "missing start"),
            Self::InvalidStart(_) => write!(f, "invalid start"),
            Self::MissingEnd => write!(f, "missing end"),
            Self::InvalidEnd(_) => write!(f, "invalid end"),
            Self::MissingScore => write!(f, "missing score"),
            Self::InvalidScore(_) => write!(f, "invalid score"),
            Self::MissingStrand => write!(f, "missing strand"),
            Self::InvalidStrand(_) => write!(f, "invalid strand"),
            Self::MissingFrame => write!(f, "missing frame"),
            Self::InvalidFrame(_) => write!(f, "invalid frame"),
            Self::MissingAttributes => write!(f, "missing attributes"),
            Self::InvalidAttributes(_) => write!(f, "invalid attributes"),
        }
    }
}

impl FromStr for RecordBuf {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        const FIELD_DELIMITER: char = '\t';
        const MAX_FIELDS: usize = 9;

        let mut fields = s.trim_end().splitn(MAX_FIELDS, FIELD_DELIMITER);

        let reference_sequence_name = fields
            .next()
            .map(|s| s.into())
            .ok_or(ParseError::MissingReferenceSequenceName)?;

        let source = fields
            .next()
            .map(|s| s.into())
            .ok_or(ParseError::MissingSource)?;

        let ty = fields
            .next()
            .map(|s| s.into())
            .ok_or(ParseError::MissingType)?;

        let start = fields
            .next()
            .ok_or(ParseError::MissingStart)
            .and_then(|s| s.parse().map_err(ParseError::InvalidStart))?;

        let end = fields
            .next()
            .ok_or(ParseError::MissingEnd)
            .and_then(|s| s.parse().map_err(ParseError::InvalidEnd))?;

        let score = fields
            .next()
            .ok_or(ParseError::MissingScore)
            .and_then(parse_score)?;

        let strand = fields
            .next()
            .ok_or(ParseError::MissingStrand)
            .and_then(parse_strand)?;

        let frame = fields
            .next()
            .ok_or(ParseError::MissingFrame)
            .and_then(parse_frame)?;

        let attributes = fields
            .next()
            .ok_or(ParseError::MissingAttributes)
            .and_then(parse_attributes)?;

        Ok(Self {
            reference_sequence_name,
            source,
            ty,
            start,
            end,
            score,
            strand,
            frame,
            attributes,
        })
    }
}

fn parse_score(s: &str) -> Result<Option<f32>, ParseError> {
    if s == MISSING_FIELD {
        Ok(None)
    } else {
        s.parse().map(Some).map_err(ParseError::InvalidScore)
    }
}

fn parse_strand(s: &str) -> Result<Option<Strand>, ParseError> {
    if s == MISSING_FIELD {
        Ok(None)
    } else {
        s.parse().map(Some).map_err(ParseError::InvalidStrand)
    }
}

fn parse_frame(s: &str) -> Result<Option<Frame>, ParseError> {
    if s == MISSING_FIELD {
        Ok(None)
    } else {
        s.parse().map(Some).map_err(ParseError::InvalidFrame)
    }
}

fn parse_attributes(s: &str) -> Result<Attributes, ParseError> {
    s.parse().map_err(ParseError::InvalidAttributes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), noodles_core::position::TryFromIntError> {
        let expected = RecordBuf {
            reference_sequence_name: String::from("sq0"),
            source: String::from("NOODLES"),
            ty: String::from("gene"),
            start: Position::try_from(8)?,
            end: Position::try_from(13)?,
            score: None,
            strand: Some(Strand::Forward),
            frame: None,
            attributes: [
                (String::from("gene_id"), String::from("g0")),
                (String::from("transcript_id"), String::from("t0")),
            ]
            .into_iter()
            .collect(),
        };

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"g0\"; transcript_id \"t0\";";
        assert_eq!(s.parse(), Ok(expected.clone()));

        let s = "sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"g0\"; transcript_id \"t0\"; ";
        assert_eq!(s.parse(), Ok(expected));

        Ok(())
    }
}
