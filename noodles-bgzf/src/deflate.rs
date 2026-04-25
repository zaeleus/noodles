use std::io;

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
    use zlib_rs::{Inflate, InflateFlush, Status};

    const HAS_ZLIB_HEADER: bool = false;
    const WINDOW_BITS: u8 = 15;

    let mut decoder = Inflate::new(HAS_ZLIB_HEADER, WINDOW_BITS);

    let status = decoder
        .decompress(src, dst, InflateFlush::Finish)
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidData))?;

    if status == Status::StreamEnd {
        Ok(())
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidData))
    }
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

    Ok(crc32(src))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn encode(src: &[u8], compression_level: i32, dst: &mut Vec<u8>) -> io::Result<u32> {
    use zlib_rs::{Deflate, DeflateFlush, Status, compress_bound};

    const HAS_ZLIB_HEADER: bool = false;
    const WINDOW_BITS: u8 = 15;

    let max_len = compress_bound(src.len());
    dst.resize(max_len, 0);

    let mut encoder = Deflate::new(compression_level, HAS_ZLIB_HEADER, WINDOW_BITS);

    let status = encoder
        .compress(src, dst, DeflateFlush::Finish)
        .map_err(|_| io::Error::from(io::ErrorKind::InvalidInput))?;

    if status == Status::StreamEnd {
        dst.truncate(encoder.total_out() as usize);
        Ok(crc32(src))
    } else {
        Err(io::Error::from(io::ErrorKind::InvalidInput))
    }
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn crc32(src: &[u8]) -> u32 {
    const START: u32 = 0;

    zlib_rs::crc32::crc32(START, src)
}

#[cfg(feature = "libdeflate")]
pub(crate) fn crc32(src: &[u8]) -> u32 {
    let mut crc = libdeflater::Crc::new();
    crc.update(src);
    crc.sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc32() {
        assert_eq!(crc32(b""), 0x00000000);
        assert_eq!(crc32(b"noodles"), 0x802a58a1);
    }
}
