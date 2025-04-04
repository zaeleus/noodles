use std::io;

use bstr::BString;
use noodles_gff::feature::RecordBuf;

use super::Line;

/// A GTF line.
#[derive(Clone, Debug, PartialEq)]
pub enum LineBuf {
    /// A comment (`#`).
    Comment(BString),
    /// A record.
    Record(RecordBuf),
}

impl TryFrom<Line> for LineBuf {
    type Error = io::Error;

    fn try_from(line: Line) -> Result<Self, Self::Error> {
        if let Some(s) = line.as_comment() {
            Ok(Self::Comment(s.into()))
        } else if let Some(result) = line.as_record() {
            result
                .and_then(|record| RecordBuf::try_from_feature_record(&record))
                .map(Self::Record)
        } else {
            unreachable!()
        }
    }
}
