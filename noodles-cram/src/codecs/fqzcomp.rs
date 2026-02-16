mod decode;
mod encode;
mod models;
mod parameters;

use self::models::Models;
pub use self::{decode::decode, encode::encode};

#[cfg(test)]
mod tests {
    use std::io;

    use super::*;

    #[test]
    fn test_self() -> io::Result<()> {
        fn t(data: &[Vec<u8>]) -> io::Result<()> {
            let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
            let src: Vec<_> = data.iter().flatten().copied().collect();

            let compressed_data = encode(&records, &src)?;
            let uncompressed_data = decode(&compressed_data)?;

            assert_eq!(uncompressed_data, src);

            Ok(())
        }

        t(&[
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0],
        ])?;

        t(&[
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        ])?;

        // HAVE_QMAP: 3 distinct values
        t(&[
            vec![5, 10, 15, 5, 10, 15, 5, 10, 15, 5],
            vec![10, 15, 5, 10, 15, 5, 10, 15, 5, 10],
            vec![15, 5, 10, 15, 5],
        ])?;

        // HAVE_QMAP: single value
        t(&[vec![42, 42, 42, 42, 42], vec![42, 42, 42, 42, 42]])?;

        // No QMAP: 17 distinct values
        t(&[
            (0u8..17).collect(),
            (0u8..17).collect(),
            (0u8..17).collect(),
        ])?;

        Ok(())
    }

    #[test]
    fn test_self_multi_param() -> io::Result<()> {
        // MULTI_PARAM: 15 short + 15 long reads
        let mut data: Vec<Vec<u8>> = Vec::new();
        for i in 0..15 {
            data.push(vec![(i % 4) as u8; 5]);
        }
        for i in 0..15 {
            data.push(vec![(i % 6) as u8; 20]);
        }

        let records: Vec<_> = data.iter().map(|scores| (scores.len(), false)).collect();
        let src: Vec<_> = data.iter().flatten().copied().collect();

        let compressed_data = encode(&records, &src)?;
        let uncompressed_data = decode(&compressed_data)?;

        assert_eq!(uncompressed_data, src);

        Ok(())
    }

    #[test]
    fn test_self_multi_param_with_reverse() -> io::Result<()> {
        // MULTI_PARAM + DO_REV
        let mut data: Vec<Vec<u8>> = Vec::new();
        let mut records = Vec::new();
        for i in 0..15 {
            data.push((0..5).map(|j| ((i * 3 + j) % 10) as u8).collect());
            records.push((5, i % 3 == 1));
        }
        for i in 0..15 {
            data.push((0..20).map(|j| ((i * 7 + j) % 12) as u8).collect());
            records.push((20, i % 4 == 0));
        }

        let src: Vec<_> = data.iter().flatten().copied().collect();

        let compressed_data = encode(&records, &src)?;
        let uncompressed_data = decode(&compressed_data)?;

        assert_eq!(uncompressed_data, src);

        Ok(())
    }
}
