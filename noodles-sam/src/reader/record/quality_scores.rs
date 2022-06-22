use std::{io, str};

use crate::record::QualityScores;

pub(crate) fn parse_quality_scores(src: &[u8]) -> io::Result<QualityScores> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(QualityScores::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}
