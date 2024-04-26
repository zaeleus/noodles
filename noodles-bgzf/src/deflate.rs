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
    use std::io::Read;

    use flate2::bufread::DeflateDecoder;

    let mut decoder = DeflateDecoder::new(src);
    decoder.read_exact(dst)
}
