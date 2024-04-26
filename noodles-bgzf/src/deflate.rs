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
) -> io::Result<(Vec<u8>, u32, u32)> {
    use libdeflater::Compressor;

    let mut encoder = Compressor::new(compression_level);

    let max_len = encoder.deflate_compress_bound(src.len());
    let mut dst = vec![0; max_len];

    let len = encoder
        .deflate_compress(src, &mut dst)
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;

    dst.resize(len, 0);

    let mut crc = Crc::new();
    crc.update(src);

    Ok((dst, crc.sum(), crc.amount()))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn encode(
    src: &[u8],
    compression_level: flate2::Compression,
) -> io::Result<(Vec<u8>, u32, u32)> {
    use flate2::write::DeflateEncoder;

    let mut encoder = DeflateEncoder::new(Vec::new(), compression_level);
    encoder.write_all(src)?;
    let dst = encoder.finish()?;

    let mut crc = Crc::new();
    crc.update(src);

    Ok((dst, crc.sum(), crc.amount()))
}
