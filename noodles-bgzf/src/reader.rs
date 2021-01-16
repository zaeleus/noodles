use std::{
    cmp,
    convert::TryFrom,
    io::{self, BufRead, Read, Seek, SeekFrom},
};

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::DeflateDecoder;

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
/// reader.read_to_end(&mut data);
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
        self.position = cpos;

        read_block(&mut self.inner, &mut self.cdata, &mut self.block)?;

        self.block.set_upos(upos as u32);

        Ok(pos)
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Same idea as in io::BufReader::Read.
        // TODO: Better handling of reads larger than an entire block?
        let bytes_read = {
            let mut remaining = self.fill_buf()?;
            remaining.read(buf)?
        };

        self.consume(bytes_read);

        Ok(bytes_read)
    }
}

impl<R> BufRead for Reader<R>
where
    R: Read,
{
    fn consume(&mut self, amt: usize) {
        // Ensure addition of amt to current upos does not does not overflow
        let upos = match u32::try_from(amt) {
            Ok(n) => match self.block.upos().checked_add(n) {
                Some(m) => m,
                None => u32::MAX,
            },
            Err(_) => u32::MAX,
        };

        let upos = cmp::min(self.block.ulen(), upos);

        self.block.set_upos(upos)
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.block.is_eof() {
            match read_block(&mut self.inner, &mut self.cdata, &mut self.block) {
                Ok(bs) => {
                    self.block.set_cpos(self.position);
                    self.position += bs as u64;
                }
                Err(e) => return Err(e),
            }
        }

        Ok(self.block.fill_buf())
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

    if reader.read_exact(&mut trailer).is_err() {
        // Block must contain valid trailer
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF trailer",
        ));
    }

    let r#isize = &trailer[4..8];

    Ok(LittleEndian::read_u32(r#isize))
}

fn inflate_data<R>(reader: R, writer: &mut Vec<u8>) -> io::Result<usize>
where
    R: Read,
{
    let mut decoder = DeflateDecoder::new(reader);
    decoder.read_to_end(writer)
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

    let cdata_len = clen - BGZF_HEADER_SIZE - gz::TRAILER_SIZE;
    cdata.resize(cdata_len, Default::default());
    reader.read_exact(cdata)?;

    let ulen = read_trailer(reader)?;

    block.set_clen(clen as u64);
    block.set_upos(0);

    let udata = block.data_mut();
    udata.clear();

    inflate_data(&cdata[..], udata)?;

    if udata.len() != ulen as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BGZF block length not equal to isize",
        ));
    }

    Ok(clen)
}

#[cfg(test)]
mod tests {
    use crate::writer::BGZF_EOF;

    use super::*;

    #[test]
    fn test_read_header() -> io::Result<()> {
        let mut reader = &BGZF_EOF[..];
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
        let mut reader = &BGZF_EOF[..];
        let mut cdata = Vec::new();
        let mut block = Block::default();

        let block_size = read_block(&mut reader, &mut cdata, &mut block)?;
        assert_eq!(block_size, BGZF_EOF.len());

        Ok(())
    }
}
