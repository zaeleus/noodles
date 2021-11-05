use std::{
    cmp,
    io::{self, BufRead, Read, Seek, SeekFrom},
};

use byteorder::{ByteOrder, LittleEndian};

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

        let block_size = read_block(&mut self.inner, &mut self.cdata, &mut self.block)?;
        self.position = cpos + (block_size as u64);

        self.block.set_cpos(cpos);
        self.block.set_upos(u32::from(upos));

        Ok(pos)
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
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
        let upos = cmp::min(self.block.ulen(), self.block.upos() + amt as u32);
        self.block.set_upos(upos);
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.block.is_eof() {
            let block_size = read_block(&mut self.inner, &mut self.cdata, &mut self.block)?;
            self.block.set_cpos(self.position);
            self.position += block_size as u64;
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

    let bsize = LittleEndian::read_u16(&header[16..]);

    // Add 1 because BSIZE is "total Block SIZE minus 1".
    Ok(u32::from(bsize) + 1)
}

/// Reads a BGZF block trailer.
///
/// The position of the stream is expected to be at the start of the block trailer, i.e., 8 bytes
/// from the end of the block.
///
/// This returns the length of the uncompressed data (`ISIZE`).
fn read_trailer<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut trailer = [0; gz::TRAILER_SIZE];
    reader.read_exact(&mut trailer)?;
    let r#isize = LittleEndian::read_u32(&trailer[4..]);
    Ok(r#isize)
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

fn read_block<R>(reader: &mut R, cdata: &mut Vec<u8>, block: &mut Block) -> io::Result<usize>
where
    R: Read,
{
    let clen = match read_header(reader) {
        Ok(0) => return Ok(0),
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
    cdata.resize(cdata_len, Default::default());
    reader.read_exact(cdata)?;

    let ulen = read_trailer(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
    })?;

    block.set_clen(clen as u64);
    block.set_upos(0);

    let udata = block.data_mut();
    udata.resize(ulen, Default::default());
    inflate_data(cdata, udata)?;

    Ok(clen)
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
    fn test_read_trailer() -> io::Result<()> {
        let (_, mut reader) = BGZF_EOF.split_at(BGZF_EOF.len() - gz::TRAILER_SIZE);
        let r#isize = read_trailer(&mut reader)?;
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
