use std::{io, mem};

use bytes::Buf;
use noodles_sam::record::QualityScores;

pub fn get_quality_scores<B>(
    src: &mut B,
    quality_scores: &mut QualityScores,
    l_seq: usize,
) -> io::Result<()>
where
    B: Buf,
{
    if l_seq == 0 {
        quality_scores.clear();
        return Ok(());
    }

    if src.remaining() < l_seq {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    if is_missing_quality_scores(src.take(l_seq).chunk()) {
        quality_scores.clear();
        src.advance(l_seq);
    } else {
        let scores = Vec::from(mem::take(quality_scores));

        let mut buf: Vec<_> = scores.into_iter().map(u8::from).collect();
        buf.resize(l_seq, 0);
        src.copy_to_slice(&mut buf);

        *quality_scores = QualityScores::try_from(buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    Ok(())
}

fn is_missing_quality_scores(src: &[u8]) -> bool {
    const MISSING: u8 = 0xff;

    src.iter().all(|&b| b == MISSING)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(mut src: &[u8], expected: &QualityScores) -> io::Result<()> {
            let mut actual = QualityScores::default();
            get_quality_scores(&mut src, &mut actual, expected.len())?;
            assert_eq!(&actual, expected);
            Ok(())
        }

        t(&[], &QualityScores::default())?;
        t(&[0x2d, 0x23, 0x2b, 0x32], &"NDLS".parse()?)?;

        Ok(())
    }

    #[test]
    fn test_get_quality_scores_with_sequence_and_no_quality_scores() -> io::Result<()> {
        let data = [0xff, 0xff, 0xff, 0xff];
        let mut buf = &data[..];

        let mut quality_scores = QualityScores::default();
        get_quality_scores(&mut buf, &mut quality_scores, 4)?;

        assert!(quality_scores.is_empty());

        Ok(())
    }
}
