use std::io;

use flate2::Crc;

#[cfg(feature = "libdeflate")]
pub(crate) fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use libdeflater::Decompressor;

    let mut decoder = Decompressor::new();

    decoder
        .deflate_decompress(src, dst)
        .map(|_| ())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn decode(src: &[u8], dst: &mut [u8]) -> io::Result<()> {
    use std::io::Read;

    use flate2::bufread::DeflateDecoder;

    let mut decoder = DeflateDecoder::new(src);
    decoder.read_exact(dst)
}

#[cfg(feature = "libdeflate")]
pub(crate) fn encode(
    src: &[u8],
    compression_level: libdeflater::CompressionLvl,
    dst: &mut Vec<u8>,
) -> io::Result<u32> {
    use libdeflater::Compressor;

    let mut encoder = Compressor::new(compression_level);

    let max_len = encoder.deflate_compress_bound(src.len());
    dst.resize(max_len, 0);

    let len = encoder
        .deflate_compress(src, dst)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    dst.truncate(len);

    let mut crc = Crc::new();
    crc.update(src);

    Ok(crc.sum())
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn encode(
    src: &[u8],
    compression_level: flate2::Compression,
    dst: &mut Vec<u8>,
) -> io::Result<u32> {
    use std::io::Write;

    use flate2::write::DeflateEncoder;

    dst.clear();

    let mut encoder = DeflateEncoder::new(dst, compression_level);
    encoder.write_all(src)?;
    encoder.finish()?;

    let mut crc = Crc::new();
    crc.update(src);

    Ok(crc.sum())
}
