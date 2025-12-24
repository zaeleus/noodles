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
            let lens: Vec<_> = data.iter().map(|scores| scores.len()).collect();
            let src: Vec<_> = data.iter().flatten().copied().collect();

            let compressed_data = encode(&lens, &src)?;
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

        Ok(())
    }
}
