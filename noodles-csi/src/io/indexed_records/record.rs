mod position;

use std::{error, fmt, ops::Range};

use noodles_core::Position;

use self::position::parse_start_position;
use crate::{
    binning_index::index::header::format::coordinate_system::CoordinateSystem, io::IndexedRecord,
};

/// An indexed record.
pub struct Record {
    buf: String,
    reference_sequence_name_bounds: Range<usize>,
    start_position: Position,
    end_position: Position,
}

impl IndexedRecord for Record {
    fn indexed_reference_sequence_name(&self) -> &str {
        &self.buf[self.reference_sequence_name_bounds.clone()]
    }

    fn indexed_start_position(&self) -> Position {
        self.start_position
    }

    fn indexed_end_position(&self) -> Position {
        self.end_position
    }
}

impl AsRef<str> for Record {
    fn as_ref(&self) -> &str {
        &self.buf
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start position is missing.
    MissingStartPosition,
    /// The start position is invalid.
    InvalidStartPosition(position::ParseError),
    /// The end position is missing.
    MissingEndPosition,
    /// The end position is invalid.
    InvalidEndPosition(position::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStartPosition(e) => Some(e),
            Self::InvalidEndPosition(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => write!(f, "missing reference sequence name"),
            Self::MissingStartPosition => write!(f, "missing start position"),
            Self::InvalidStartPosition(_) => write!(f, "invalid start position"),
            Self::MissingEndPosition => write!(f, "missing end position"),
            Self::InvalidEndPosition(_) => write!(f, "invalid end position"),
        }
    }
}

pub(crate) fn parse_record(
    s: String,
    reference_sequence_name_index: usize,
    start_position_index: usize,
    end_position_index: Option<usize>,
    coordinate_system: CoordinateSystem,
) -> Result<Record, ParseError> {
    const DELIMITER: char = '\t';

    let fields: Vec<_> = s.split(DELIMITER).collect();

    let reference_sequence_name_bounds =
        calculate_reference_sequence_name_bounds(&fields, reference_sequence_name_index)?;

    let raw_start = fields
        .get(start_position_index)
        .ok_or(ParseError::MissingStartPosition)?;

    let start_position = parse_start_position(raw_start, coordinate_system)
        .map_err(ParseError::InvalidStartPosition)?;

    let end_position = if let Some(i) = end_position_index {
        fields
            .get(i)
            .ok_or(ParseError::MissingEndPosition)
            .and_then(|s| {
                s.parse()
                    .map_err(position::ParseError::Parse)
                    .map_err(ParseError::InvalidEndPosition)
            })?
    } else {
        // _The Tabix index file format_: "Field `col_beg` may equal `col_end`, and in this case,
        // the end of a region is `end=beg+1`."
        start_position
            .checked_add(1)
            .expect("attempt to add with overflow")
    };

    Ok(Record {
        buf: s,
        reference_sequence_name_bounds,
        start_position,
        end_position,
    })
}

fn calculate_reference_sequence_name_bounds(
    fields: &[&str],
    i: usize,
) -> Result<Range<usize>, ParseError> {
    let name = fields
        .get(i)
        .ok_or(ParseError::MissingReferenceSequenceName)?;

    let mut start = 0;

    for s in fields.iter().take(i) {
        start += s.len() + 1;
    }

    let end = start + name.len();

    Ok(start..end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_record() -> Result<(), ParseError> {
        let record = parse_record("sq0\t7\t13".into(), 0, 1, Some(2), CoordinateSystem::Bed)?;
        assert_eq!(record.reference_sequence_name_bounds, 0..3);
        assert_eq!(record.start_position, const { Position::new(8).unwrap() });
        assert_eq!(record.end_position, const { Position::new(13).unwrap() });

        let record = parse_record("sq0\t8\t13".into(), 0, 1, Some(2), CoordinateSystem::Gff)?;
        assert_eq!(record.reference_sequence_name_bounds, 0..3);
        assert_eq!(record.start_position, const { Position::new(8).unwrap() });
        assert_eq!(record.end_position, const { Position::new(13).unwrap() });

        let record = parse_record("sq0\t8".into(), 0, 1, None, CoordinateSystem::Gff)?;
        assert_eq!(record.reference_sequence_name_bounds, 0..3);
        assert_eq!(record.start_position, const { Position::new(8).unwrap() });
        // FIXME: `end_position` = 8.
        assert_eq!(record.end_position, const { Position::new(9).unwrap() });

        // FIXME: `ParseError::MissingReferenceSequenceName`.
        assert!(matches!(
            parse_record("".into(), 0, 1, Some(2), CoordinateSystem::Gff),
            Err(ParseError::MissingStartPosition)
        ));

        assert!(matches!(
            parse_record("sq0".into(), 0, 1, Some(2), CoordinateSystem::Gff),
            Err(ParseError::MissingStartPosition)
        ));

        assert!(matches!(
            parse_record("sq0\tn".into(), 0, 1, Some(2), CoordinateSystem::Gff),
            Err(ParseError::InvalidStartPosition(_))
        ));

        assert!(matches!(
            parse_record("sq0\t8".into(), 0, 1, Some(2), CoordinateSystem::Gff),
            Err(ParseError::MissingEndPosition)
        ));

        assert!(matches!(
            parse_record("sq0\t8\tn".into(), 0, 1, Some(2), CoordinateSystem::Gff),
            Err(ParseError::InvalidEndPosition(_))
        ));

        Ok(())
    }
}
