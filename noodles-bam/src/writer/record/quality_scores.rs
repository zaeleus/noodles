use bytes::BufMut;
use noodles_sam::alignment::record::AlignmentQualityScores;

pub fn put_quality_scores<B, S>(dst: &mut B, quality_scores: &S)
where
    B: BufMut,
    S: AlignmentQualityScores + ?Sized,
{
    for score in quality_scores.scores() {
        dst.put_u8(u8::from(score));
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam as sam;

    use super::*;

    #[test]
    fn test_put_quality_scores() -> Result<(), sam::alignment::record::quality_scores::ParseError> {
        use sam::alignment::record::QualityScores;

        fn t(buf: &mut Vec<u8>, quality_scores: &QualityScores, expected: &[u8]) {
            buf.clear();
            put_quality_scores(buf, quality_scores);
            assert_eq!(buf, expected);
        }

        let mut buf = Vec::new();

        t(&mut buf, &QualityScores::default(), &[]);
        t(&mut buf, &"NDLS".parse()?, &[0x2d, 0x23, 0x2b, 0x32]);

        Ok(())
    }
}
