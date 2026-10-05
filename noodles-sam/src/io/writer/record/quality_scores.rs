use std::io::{self, Write};

use super::MISSING;
use crate::alignment::record::{QualityScores, QualityScoresRef};

const OFFSET: u8 = b'!';

pub(super) fn write_quality_scores<W>(
    writer: &mut W,
    base_count: usize,
    quality_scores: QualityScoresRef<'_>,
) -> io::Result<()>
where
    W: Write,
{
    if quality_scores.is_empty() {
        writer.write_all(&[MISSING])?;
    } else if quality_scores.len() == base_count {
        match quality_scores {
            QualityScoresRef::Raw(s) => write_raw_quality_scores(writer, s)?,
            QualityScoresRef::Offset(s, offset) => write_offset_quality_scores(writer, s, offset)?,
            QualityScoresRef::QualityScores(s) => write_generic_quality_scores(writer, s)?,
        }
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

fn write_raw_quality_scores<W>(writer: &mut W, quality_scores: &[u8]) -> io::Result<()>
where
    W: Write,
{
    const CHUNK_SIZE: usize = 256;

    let mut buf = [0; CHUNK_SIZE];

    if quality_scores.iter().all(|&n| is_valid_score(n)) {
        for chunk in quality_scores.chunks(CHUNK_SIZE) {
            for (&s, d) in chunk.iter().zip(&mut buf) {
                *d = encode(s);
            }

            let len = chunk.len();
            writer.write_all(&buf[..len])?;
        }

        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
}

fn write_offset_quality_scores<W>(
    writer: &mut W,
    quality_scores: &[u8],
    offset: u8,
) -> io::Result<()>
where
    W: Write,
{
    match offset {
        0 => write_raw_quality_scores(writer, quality_scores)?,
        OFFSET => {
            if quality_scores.iter().all(|&n| is_valid_encoded_score(n)) {
                writer.write_all(quality_scores)?;
            } else {
                return Err(io::Error::from(io::ErrorKind::InvalidInput));
            }
        }
        _ => {
            for n in quality_scores {
                let m = n
                    .checked_sub(offset)
                    .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidInput))?;

                if is_valid_score(m) {
                    writer.write_all(&[encode(m)])?;
                } else {
                    return Err(io::Error::from(io::ErrorKind::InvalidInput));
                }
            }
        }
    }

    Ok(())
}

fn write_generic_quality_scores<W, S>(writer: &mut W, quality_scores: S) -> io::Result<()>
where
    W: Write,
    S: QualityScores,
{
    for result in quality_scores.iter() {
        let n = result?;

        if is_valid_score(n) {
            let m = encode(n);
            writer.write_all(&[m])?;
        } else {
            return Err(io::Error::from(io::ErrorKind::InvalidInput));
        }
    }

    Ok(())
}

fn is_valid_score(score: u8) -> bool {
    const MAX_SCORE: u8 = b'~' - OFFSET;
    score <= MAX_SCORE
}

fn is_valid_encoded_score(n: u8) -> bool {
    const MIN: u8 = OFFSET;
    const MAX: u8 = b'~';

    (MIN..=MAX).contains(&n)
}

fn encode(n: u8) -> u8 {
    n + OFFSET
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::record_buf::QualityScores as QualityScoresBuf;

    #[test]
    fn test_write_quality_scores() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            buf: &mut Vec<u8>,
            base_count: usize,
            quality_scores: &QualityScoresBuf,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            let s = QualityScoresRef::QualityScores(Box::new(quality_scores));
            write_quality_scores(buf, base_count, s)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, 0, &QualityScoresBuf::default(), b"*")?;
        t(&mut buf, 4, &QualityScoresBuf::default(), b"*")?;

        let quality_scores = [45, 35, 43, 50].into_iter().collect();
        t(&mut buf, 4, &quality_scores, b"NDLS")?;

        buf.clear();
        let quality_scores = [45, 35, 43, 50].into_iter().collect();
        let s = QualityScoresRef::QualityScores(Box::new(&quality_scores));
        assert!(matches!(
            write_quality_scores(&mut buf, 3, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let quality_scores = [255].into_iter().collect();
        let s = QualityScoresRef::QualityScores(Box::new(&quality_scores));
        assert!(matches!(
            write_quality_scores(&mut buf, 1, s),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_offset_quality_scores() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_offset_quality_scores(&mut buf, &[45, 35, 43, 50], 0)?;
        assert_eq!(buf, b"NDLS");

        buf.clear();
        write_offset_quality_scores(&mut buf, b"NDLS", b'!')?;
        assert_eq!(buf, b"NDLS");

        buf.clear();
        write_offset_quality_scores(&mut buf, &[46, 36, 44, 51], 1)?;
        assert_eq!(buf, b"NDLS");

        buf.clear();
        let quality_scores = [255];
        assert!(matches!(
            write_offset_quality_scores(&mut buf, &quality_scores, OFFSET),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_raw_quality_scores() -> io::Result<()> {
        fn t(buf: &mut Vec<u8>, src: &[u8], expected: &[u8]) -> io::Result<()> {
            buf.clear();
            write_raw_quality_scores(buf, src)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &[45, 35, 43, 50], b"NDLS")?;

        let src = [45, 35, 43, 50].repeat(96);
        let expected = b"NDLS".repeat(96);
        t(&mut buf, &src, &expected)?;

        buf.clear();
        let quality_scores = [255];
        assert!(matches!(
            write_raw_quality_scores(&mut buf, &quality_scores),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(0), b'!');
        assert_eq!(encode(93), b'~');
    }
}
