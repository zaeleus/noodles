use crate::record::{alternate_bases, AlternateBases};

pub(super) fn parse_alternate_bases(
    s: &str,
) -> Result<AlternateBases, alternate_bases::ParseError> {
    s.parse()
}
