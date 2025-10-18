use std::io::{self, Read, Write};

use liblzma::{bufread::XzDecoder, write::XzEncoder};

pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    let mut decoder = XzDecoder::new(src);
    decoder.read_exact(dst)
}

pub fn encode(compression_level: u32, src: &[u8]) -> io::Result<Vec<u8>> {
    let mut encoder = XzEncoder::new(Vec::new(), compression_level);
    encoder.write_all(src)?;
    encoder.finish()
}
