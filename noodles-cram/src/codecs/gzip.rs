use std::io;

use flate2::Compression;

#[cfg(feature = "libdeflate")]
pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use libdeflater::Decompressor;

    let mut decoder = Decompressor::new();

    decoder
        .gzip_decompress(src, dst)
        .map(|_| ())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(not(feature = "libdeflate"))]
pub fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use std::io::Read;

    use flate2::bufread::GzDecoder;

    let mut decoder = GzDecoder::new(src);
    decoder.read_exact(dst)
}

#[cfg(feature = "libdeflate")]
pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    use libdeflater::{CompressionLvl, Compressor};

    let lvl = i32::try_from(compression_level.level())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
        .and_then(|n| {
            CompressionLvl::new(n).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid libdeflate compression level",
                )
            })
        })?;

    let mut encoder = Compressor::new(lvl);

    let max_len = encoder.gzip_compress_bound(src.len());
    let mut dst = vec![0; max_len];

    let len = encoder
        .gzip_compress(src, &mut dst)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    dst.resize(len, 0);

    Ok(dst)
}

#[cfg(not(feature = "libdeflate"))]
pub fn encode(compression_level: Compression, src: &[u8]) -> io::Result<Vec<u8>> {
    use std::io::Write;

    use flate2::write::GzEncoder;

    let mut encoder = GzEncoder::new(Vec::new(), compression_level);
    encoder.write_all(src)?;
    encoder.finish()
}
