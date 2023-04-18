use std::io;

use crate::record::MappingQuality;

pub(crate) fn parse_mapping_quality(src: &[u8]) -> io::Result<Option<MappingQuality>> {
    lexical_core::parse(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(MappingQuality::new)
}
