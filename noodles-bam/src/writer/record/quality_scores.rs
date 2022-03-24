use bytes::BufMut;
use noodles_sam as sam;

pub fn put_quality_scores<B>(dst: &mut B, quality_scores: &sam::record::QualityScores)
where
    B: BufMut,
{
    for &score in quality_scores.as_ref() {
        dst.put_u8(u8::from(score));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_quality_scores() -> Result<(), sam::record::quality_scores::ParseError> {
        fn t(buf: &mut Vec<u8>, quality_scores: &sam::record::QualityScores, expected: &[u8]) {
            buf.clear();
            put_quality_scores(buf, quality_scores);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, &sam::record::QualityScores::default(), &[]);
        t(&mut buf, &"NDLS".parse()?, &[0x2d, 0x23, 0x2b, 0x32]);

        Ok(())
    }
}
