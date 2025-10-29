use std::io::{self, Read, Write};

use lzma_rust2::{XzOptions, XzReader, XzWriter};

pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    let mut decoder = XzReader::new(src, false);
    decoder.read_exact(dst)
}

pub fn encode(compression_level: u32, src: &[u8]) -> io::Result<Vec<u8>> {
    let options = XzOptions::with_preset(compression_level);
    let mut encoder = XzWriter::new(Vec::new(), options)?;
    encoder.write_all(src)?;
    encoder.finish()
}
