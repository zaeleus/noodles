use std::io;

use bytes::BufMut;
use noodles_sam::record::QualityScores;

pub fn put_quality_scores<B>(
    dst: &mut B,
    base_count: usize,
    quality_scores: &QualityScores,
) -> io::Result<()>
where
    B: BufMut,
{
    // § 4.2.3 SEQ and QUAL encoding (2022-08-22)
    const MISSING: u8 = 255;

    if quality_scores.len() == base_count {
        for &score in quality_scores.as_ref() {
            dst.put_u8(u8::from(score));
        }
    } else if quality_scores.is_empty() {
        dst.put_bytes(MISSING, base_count);
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "sequence-quality scores length mismatch: expected {}, got {}",
                base_count,
                quality_scores.len()
            ),
        ));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_put_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            base_count: usize,
            quality_scores: &QualityScores,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            put_quality_scores(buf, base_count, quality_scores)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &QualityScores::default(), &[])?;
        t(
            &mut buf,
            4,
            &QualityScores::default(),
            &[0xff, 0xff, 0xff, 0xff],
        )?;
        t(&mut buf, 4, &"NDLS".parse()?, &[0x2d, 0x23, 0x2b, 0x32])?;

        let quality_scores = "NDLS".parse()?;
        buf.clear();
        assert!(matches!(
            put_quality_scores(&mut buf, 3, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }
}
