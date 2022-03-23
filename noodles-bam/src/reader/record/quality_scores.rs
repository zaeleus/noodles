use std::io;

use bytes::Buf;
use noodles_sam as sam;

pub(super) fn get_quality_scores<B>(
    buf: &mut B,
    quality_scores: &mut sam::record::QualityScores,
    l_seq: usize,
) -> io::Result<()>
where
    B: Buf,
{
    use sam::record::quality_scores::Score;

    if buf.remaining() < l_seq {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    quality_scores.clear();

    let qual = buf.take(l_seq);

    if !is_missing_quality_scores(qual.chunk()) {
        for &b in qual.chunk() {
            let score =
                Score::try_from(b).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            quality_scores.push(score);
        }
    }

    buf.advance(l_seq);

    Ok(())
}

fn is_missing_quality_scores(buf: &[u8]) -> bool {
    use crate::writer::alignment_record::NULL_QUALITY_SCORE;

    buf.iter().all(|&b| b == NULL_QUALITY_SCORE)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_quality_scores_with_sequence_and_no_quality_scores() -> io::Result<()> {
        let data = [0xff, 0xff, 0xff, 0xff];
        let mut buf = &data[..];

        let mut quality_scores = sam::record::QualityScores::default();
        get_quality_scores(&mut buf, &mut quality_scores, 4)?;

        assert!(quality_scores.is_empty());

        Ok(())
    }
}
