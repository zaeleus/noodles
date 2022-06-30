use std::{
    cmp,
    io::{self, BufRead, Read, Seek, SeekFrom},
};

use byteorder::{ByteOrder, LittleEndian, ReadBytesExt};
use flate2::Crc;

use super::{gz, Block, VirtualPosition, BGZF_HEADER_SIZE};

/// A BGZF reader.
///
/// Due to the static structure of a BGZF block, gzip headers are mostly discarded. CRC32
/// validation is also disabled when decompressing data.
///
/// The reader implements both [`std::io::Read`] and [`std::io::BufRead`], consuming compressed
/// data and emitting uncompressed data. It is internally buffered by a single block, and to
/// correctly track (virtual) positions, the reader _cannot_ be double buffered (e.g., using
/// [`std::io::BufReader`]).
///
/// # Examples
///
/// ```no_run
/// # use std::{fs::File, io::{self, Read}};
/// use noodles_bgzf as bgzf;
/// let mut reader = File::open("data.gz").map(bgzf::Reader::new)?;
/// let mut data = Vec::new();
/// reader.read_to_end(&mut data)?;
/// # Ok::<(), io::Error>(())
/// ```
pub struct Reader<R> {
    inner: R,
    position: u64,
    cdata: Vec<u8>,
    block: Block,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            position: 0,
            cdata: Vec::new(),
            block: Block::default(),
        }
    }

    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let mut reader = bgzf::Reader::new(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Unwraps and returns the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// Returns the current position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// assert_eq!(reader.position(), 0);
    /// ```
    pub fn position(&self) -> u64 {
        self.position
    }

    /// Returns the current virtual position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// assert_eq!(reader.virtual_position(), bgzf::VirtualPosition::from(0));
    /// ```
    pub fn virtual_position(&self) -> VirtualPosition {
        self.block.virtual_position()
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the stream to the given virtual position.
    ///
    /// The underlying stream's cursor is first moved the the compressed position. A block is read,
    /// decompressed, and has its own cursor moved to the uncompressed position.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io::{self, Cursor};
    /// use noodles_bgzf as bgzf;
    /// let mut reader = bgzf::Reader::new(Cursor::new(Vec::new()));
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let (cpos, upos) = pos.into();

        self.inner.seek(SeekFrom::Start(cpos))?;

        read_block(&mut self.inner, &mut self.cdata, &mut self.block)?;
        self.position = cpos + self.block.clen();

        self.block.set_cpos(cpos);
        self.block.set_upos(usize::from(upos));

        Ok(pos)
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        use super::block;

        // If a new block is about to be read and the given buffer is guaranteed to be larger than
        // next block, reading to the block buffer can be skipped. The uncompressed data is read
        // directly to the given buffer to avoid double copying.
        if self.block.is_eof() && buf.len() >= block::MAX_UNCOMPRESSED_DATA_LENGTH {
            read_block_into(&mut self.inner, &mut self.cdata, &mut self.block, buf)?;
            self.block.set_cpos(self.position);
            self.position += self.block.clen();
            return Ok(self.block.ulen());
        }

        let bytes_read = {
            let mut remaining = self.fill_buf()?;
            remaining.read(buf)?
        };

        self.consume(bytes_read);

        Ok(bytes_read)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        let remaining = self.block.buffer();

        if buf.len() <= remaining.len() {
            buf.copy_from_slice(&remaining[..buf.len()]);
            self.consume(buf.len());
            Ok(())
        } else {
            default_read_exact(self, buf)
        }
    }
}

impl<R> BufRead for Reader<R>
where
    R: Read,
{
    fn consume(&mut self, mut amt: usize) {
        amt = cmp::min(amt, crate::block::MAX_UNCOMPRESSED_DATA_LENGTH);
        let upos = cmp::min(self.block.ulen(), self.block.upos() + amt);
        self.block.set_upos(upos);
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.block.is_eof() {
            read_block(&mut self.inner, &mut self.cdata, &mut self.block)?;
            self.block.set_cpos(self.position);
            self.position += self.block.clen();
        }

        Ok(self.block.buffer())
    }
}

/// Reads a BGZF block header.
///
/// The position of the stream is expected to be at the start of a block.
///
/// If successful, the block size (`BSIZE` + 1) is returned. If a block size of 0 is returned, the
/// stream reached EOF.
fn read_header<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut header = [0; BGZF_HEADER_SIZE];

    match reader.read_exact(&mut header) {
        Ok(_) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    }

    if !is_valid_header(&header) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }

    let bsize = LittleEndian::read_u16(&header[16..]);

    // Add 1 because BSIZE is "total Block SIZE minus 1".
    Ok(u32::from(bsize) + 1)
}

fn is_valid_header(header: &[u8; BGZF_HEADER_SIZE]) -> bool {
    const BGZF_CM: u8 = 0x08; // DEFLATE
    const BGZF_FLG: u8 = 0x04; // FEXTRA
    const BGZF_XLEN: u16 = 6;
    const BGZF_SI1: u8 = b'B';
    const BGZF_SI2: u8 = b'C';
    const BGZF_SLEN: u16 = 2;

    let magic_number = &header[0..2];
    let cm = header[2];
    let flg = header[3];
    let xlen = LittleEndian::read_u16(&header[10..]);
    let subfield_id_1 = header[12];
    let subfield_id_2 = header[13];
    let bc_len = LittleEndian::read_u16(&header[14..]);

    magic_number == gz::MAGIC_NUMBER
        && cm == BGZF_CM
        && flg == BGZF_FLG
        && xlen == BGZF_XLEN
        && subfield_id_1 == BGZF_SI1
        && subfield_id_2 == BGZF_SI2
        && bc_len == BGZF_SLEN
}

/// Reads a BGZF block trailer.
///
/// The position of the stream is expected to be at the start of the block trailer, i.e., 8 bytes
/// from the end of the block.
fn read_trailer<R>(reader: &mut R) -> io::Result<(u32, usize)>
where
    R: Read,
{
    let crc32 = reader.read_u32::<LittleEndian>()?;

    let r#isize = reader.read_u32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    Ok((crc32, r#isize))
}

#[cfg(feature = "libdeflate")]
pub(crate) fn inflate_data(reader: &[u8], writer: &mut [u8]) -> io::Result<()> {
    use libdeflater::Decompressor;

    let mut decoder = Decompressor::new();

    decoder
        .deflate_decompress(reader, writer)
        .map(|_| ())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(not(feature = "libdeflate"))]
pub(crate) fn inflate_data(reader: &[u8], writer: &mut [u8]) -> io::Result<()> {
    use flate2::bufread::DeflateDecoder;

    let mut decoder = DeflateDecoder::new(reader);
    decoder.read_exact(writer)
}

fn read_compressed_block<R>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<(usize, (u32, usize))>
where
    R: Read,
{
    let clen = match read_header(reader) {
        Ok(0) => return Ok((0, (0, 0))),
        Ok(bs) => bs as usize,
        Err(e) => return Err(e),
    };

    if clen < BGZF_HEADER_SIZE + gz::TRAILER_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected clen >= {}, got {}",
                BGZF_HEADER_SIZE + gz::TRAILER_SIZE,
                clen
            ),
        ));
    }

    let cdata_len = clen - BGZF_HEADER_SIZE - gz::TRAILER_SIZE;
    buf.resize(cdata_len, Default::default());
    reader.read_exact(buf)?;

    let trailer = read_trailer(reader)?;

    Ok((clen, trailer))
}

fn read_block<R>(reader: &mut R, cdata: &mut Vec<u8>, block: &mut Block) -> io::Result<usize>
where
    R: Read,
{
    let (clen, crc32, ulen) = match read_compressed_block(reader, cdata) {
        Ok((0, (_, 0))) => return Ok(0),
        Ok((clen, (crc32, ulen))) => (clen, crc32, ulen),
        Err(e) => return Err(e),
    };

    block.set_clen(clen as u64);
    block.set_upos(0);
    block.set_ulen(ulen);

    inflate_data(cdata, block.buffer_mut())?;

    let mut crc = Crc::new();
    crc.update(block.buffer());

    if crc.sum() == crc32 {
        Ok(clen)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}

fn read_block_into<R>(
    reader: &mut R,
    cdata: &mut Vec<u8>,
    block: &mut Block,
    buf: &mut [u8],
) -> io::Result<usize>
where
    R: Read,
{
    let (clen, crc32, ulen) = match read_compressed_block(reader, cdata) {
        Ok((0, (_, 0))) => return Ok(0),
        Ok((clen, (crc32, ulen))) => (clen, crc32, ulen),
        Err(e) => return Err(e),
    };

    block.set_clen(clen as u64);
    block.set_upos(ulen);
    block.set_ulen(ulen);

    inflate_data(cdata, &mut buf[..ulen])?;

    let mut crc = Crc::new();
    crc.update(&buf[..ulen]);

    if crc.sum() == crc32 {
        Ok(clen)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ))
    }
}

/// This is effectively the same as `std::io::default_read_exact`.
fn default_read_exact<R>(reader: &mut R, mut buf: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    while !buf.is_empty() {
        match reader.read(buf) {
            Ok(0) => break,
            Ok(n) => buf = &mut buf[n..],
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    if buf.is_empty() {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "failed to fill whole buffer",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use crate::writer::BGZF_EOF;

    use super::*;

    #[test]
    fn test_seek() -> Result<(), Box<dyn std::error::Error>> {
        #[rustfmt::skip]
        let data = [
            // block 0
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let eof = VirtualPosition::try_from((63, 0))?;

        let mut reader = Reader::new(Cursor::new(&data));

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(reader.virtual_position(), eof);

        reader.seek(VirtualPosition::try_from((0, 3))?)?;

        buf.clear();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"dles");
        assert_eq!(reader.virtual_position(), eof);

        Ok(())
    }

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = BGZF_EOF;
        let block_size = read_header(&mut reader)?;
        assert_eq!(block_size, BGZF_EOF.len() as u32);
        Ok(())
    }

    #[test]
    fn test_is_valid_header() {
        let mut src = [
            0x1f, 0x8b, // ID1, ID2
            0x08, // CM = DEFLATE
            0x04, // FLG = FEXTRA
            0x00, 0x00, 0x00, 0x00, // MTIME = 0
            0x00, // XFL = 0
            0xff, // OS = 255 (unknown)
            0x06, 0x00, // XLEN = 6
            b'B', b'C', // SI1, SI2
            0x02, 0x00, // SLEN = 2
            0x1b, 0x00, // BSIZE = 27
        ];

        assert!(is_valid_header(&src));

        src[0] = 0x00;
        assert!(!is_valid_header(&src));
    }

    #[test]
    fn test_read_trailer() -> io::Result<()> {
        let (_, mut reader) = BGZF_EOF.split_at(BGZF_EOF.len() - gz::TRAILER_SIZE);

        let (crc32, r#isize) = read_trailer(&mut reader)?;
        assert_eq!(crc32, 0);
        assert_eq!(r#isize, 0);

        Ok(())
    }

    #[test]
    fn test_read_trailer_with_invalid_block_trailer() {
        let data = [0, 0, 0, 0];
        let mut reader = &data[..];
        assert!(read_trailer(&mut reader).is_err());
    }

    #[test]
    fn test_read_block() -> io::Result<()> {
        let mut reader = BGZF_EOF;
        let mut cdata = Vec::new();
        let mut block = Block::default();

        let block_size = read_block(&mut reader, &mut cdata, &mut block)?;
        assert_eq!(block_size, BGZF_EOF.len());

        Ok(())
    }

    #[test]
    fn test_read_block_with_invalid_block_size() {
        let data = {
            let mut eof = BGZF_EOF.to_vec();
            // BSIZE = 0
            eof[16] = 0x00;
            eof[17] = 0x00;
            eof
        };

        let mut reader = &data[..];
        let mut cdata = Vec::new();
        let mut block = Block::default();

        assert!(read_block(&mut reader, &mut cdata, &mut block).is_err());
    }
}
