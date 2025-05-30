use std::io::{self, Read, Write};

use bzip2::{Compression, read::BzDecoder, write::BzEncoder};

pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    let mut decoder = BzDecoder::new(src);
    decoder.read_exact(dst)
}

pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut encoder = BzEncoder::new(Vec::new(), compression_level);
    encoder.write_all(src)?;
    encoder.finish()
}
