use std::io::{self, Read, Write};

use xz2::{bufread::XzDecoder, write::XzEncoder};

const DEFAULT_LZMA_COMPRESSION_LEVEL: u32 = 6;

pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    let mut decoder = XzDecoder::new(src);
    decoder.read_exact(dst)
}

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    let mut encoder = XzEncoder::new(Vec::new(), DEFAULT_LZMA_COMPRESSION_LEVEL);
    encoder.write_all(src)?;
    encoder.finish()
}
