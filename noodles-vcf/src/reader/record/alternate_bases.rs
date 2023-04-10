use crate::record::{alternate_bases, AlternateBases};

use super::MISSING;

pub(super) fn parse_alternate_bases(
    s: &str,
) -> Result<AlternateBases, alternate_bases::ParseError> {
    match s {
        MISSING => Ok(AlternateBases::default()),
        _ => s.parse(),
    }
}
