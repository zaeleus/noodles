//! BGZF reader.

pub(crate) mod block;
mod builder;

pub use self::builder::Builder;

use std::io::{self, BufRead, Read, Seek, SeekFrom};

use super::{gzi, Block, VirtualPosition};

/// A BGZF reader.
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
    inner: block::Inner<R>,
    position: u64,
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
        Builder::default().build_from_reader(inner)
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
        self.inner.get_ref()
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
        self.inner.get_mut()
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
        self.inner.into_inner()
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

    fn read_block(&mut self) -> io::Result<()> {
        if let Some(mut block) = self.inner.next_block()? {
            block.set_position(self.position);
            self.position += block.size();
            self.block = block;
        }

        Ok(())
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

        self.inner.get_mut().seek(SeekFrom::Start(cpos))?;
        self.position = cpos;

        self.read_block()?;

        self.block.data_mut().set_position(usize::from(upos));

        Ok(pos)
    }

    /// Seeks the stream to the given uncompressed position.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_bgzf as bgzf;
    /// let mut reader = bgzf::Reader::new(io::empty());
    /// let index = vec![(0, 0)];
    /// reader.seek_by_uncompressed_position(&index, 0)?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn seek_by_uncompressed_position(
        &mut self,
        index: &gzi::Index,
        pos: u64,
    ) -> io::Result<u64> {
        assert!(!index.is_empty());

        let i = index.iter().position(|r| r.1 > pos).unwrap_or(index.len());
        // SAFETY: `i` is > 0.
        let record = index[i - 1];

        let cpos = record.0;
        self.inner.get_mut().seek(SeekFrom::Start(cpos))?;
        self.position = cpos;

        self.read_block()?;

        let upos = usize::try_from(pos - record.1)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        self.block.data_mut().set_position(upos);

        Ok(pos)
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for Reader<R>
where
    R: Read,
{
    fn consume(&mut self, amt: usize) {
        self.block.data_mut().consume(amt)
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.block.data().has_remaining() {
            self.read_block()?;
        }

        Ok(self.block.data().as_ref())
    }
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

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_seek() -> Result<(), Box<dyn std::error::Error>> {
        #[rustfmt::skip]
        let data = [
            // block 0 (b"noodles")
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
    fn test_seek_by_uncompressed_position() -> io::Result<()> {
        #[rustfmt::skip]
        let data = [
            // block 0 (b"noodles")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x22, 0x00, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0x06, 0x00, 0xa1,
            0x58, 0x2a, 0x80, 0x07, 0x00, 0x00, 0x00,
            // block 1 (b"bgzf")
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1f, 0x00, 0x4b, 0x4a, 0xaf, 0x4a, 0x03, 0x00, 0x20, 0x68, 0xf2, 0x8c,
            0x04, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let index = vec![(0, 0), (35, 7)];

        let mut reader = Reader::new(Cursor::new(&data));

        reader.seek_by_uncompressed_position(&index, 3)?;
        let mut buf = [0; 4];
        reader.read_exact(&mut buf)?;
        assert_eq!(&buf, b"dles");

        reader.seek_by_uncompressed_position(&index, 8)?;
        let mut buf = [0; 2];
        reader.read_exact(&mut buf)?;
        assert_eq!(&buf, b"gz");

        Ok(())
    }
}
