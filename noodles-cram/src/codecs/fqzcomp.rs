mod decode;
mod encode;
mod parameter;
mod parameters;

pub use self::{decode::decode, encode::encode};

use super::aac::Model;

struct Models {
    len: Vec<Model>,
    qual: Vec<Model>,
    dup: Model,
    rev: Model,
    sel: Model,
}

impl Models {
    fn new(max_sym: u8, max_sel: u8) -> Models {
        Self {
            len: vec![Model::new(u8::MAX); 4],
            qual: vec![Model::new(max_sym); 1 << 16],
            dup: Model::new(1),
            rev: Model::new(1),
            sel: Model::new(max_sel),
        }
    }
}

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
