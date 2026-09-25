use std::{error, fmt};

use noodles_core::Position;

/// An error returned when the indexer fails to add a record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AddRecordError {
    /// The interval is invalid.
    InvalidInterval {
        /// The interval start position.
        start: Position,
        /// The interval end position.
        end: Position,
    },
    /// The end exceeds the max supported position.
    EndExceedsMaxPosition {
        /// The alignment context end position.
        end: Position,
        /// The max supported position.
        max_position: Position,
    },
    /// The reference sequence ID is out of order.
    OutOfOrderReferenceSequenceId {
        /// The given reference sequence ID.
        reference_sequence_id: usize,
        /// The current reference sequence ID.
        current_reference_sequence_id: usize,
    },
    /// The start is out of order.
    OutOfOrderStartPosition {
        /// The alignment context start position.
        start: Position,
        /// The previous start position.
        prev_start_position: Position,
    },
}

impl error::Error for AddRecordError {}

impl fmt::Display for AddRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInterval { start, end } => {
                write!(f, "invalid interval: start ({start}) > end ({end})")
            }
            Self::EndExceedsMaxPosition { end, max_position } => write!(
                f,
                "end exceeds max position: end ({end}) > max position ({max_position})"
            ),
            Self::OutOfOrderReferenceSequenceId {
                reference_sequence_id,
                current_reference_sequence_id,
            } => write!(
                f,
                "out of order reference sequence ID: reference sequence ID ({reference_sequence_id}) < current reference sequence ID ({current_reference_sequence_id})"
            ),
            Self::OutOfOrderStartPosition {
                start,
                prev_start_position,
            } => write!(
                f,
                "out of order start position: start ({start}) < previous start position ({prev_start_position})"
            ),
        }
    }
}
