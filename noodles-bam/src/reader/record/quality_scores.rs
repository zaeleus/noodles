use std::io;

use bytes::{Buf, BytesMut};

use crate::record::QualityScores;

pub fn get_quality_scores(
    src: &mut BytesMut,
    quality_scores: &mut QualityScores,
    l_seq: usize,
) -> io::Result<()> {
    use crate::reader::alignment_record::quality_scores::is_missing_quality_scores;

    if src.remaining() < l_seq {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    quality_scores.clear();

    let buf = src.split_to(l_seq);

    if !is_missing_quality_scores(&buf) {
        quality_scores.buf = buf;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_quality_scores() -> io::Result<()> {
        let data = [0x2d, 0x23, 0x2b, 0x32];
        let mut reader = BytesMut::from(&data[..]);
        let mut quality_scores = QualityScores::default();
        get_quality_scores(&mut reader, &mut quality_scores, data.len())?;
        assert_eq!(&quality_scores.buf[..], data);

        let data = [0xff, 0xff, 0xff, 0xff];
        let mut reader = BytesMut::from(&data[..]);
        get_quality_scores(&mut reader, &mut quality_scores, data.len())?;
        assert!(quality_scores.buf.is_empty());

        Ok(())
    }
}
