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
    gzi: Option<gzi::Index>,
    uncompressed_position: Option<u64>,
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

    /// Creates a BGZF reader with GZI.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let gzi = vec![];
    /// let reader = bgzf::Reader::with_gzi(&data[..], gzi);
    /// ```
    pub fn with_gzi(inner: R, gzi: gzi::Index) -> Self {
        let mut reader = Reader::new(inner);
        reader.gzi = Some(gzi);
        reader.uncompressed_position = Some(0);
        reader
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

    /// Returns the current uncompressed position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let reader = bgzf::Reader::new(&data[..]);
    /// assert_eq!(reader.uncompressed_position(), None);
    /// let gzi = vec![];
    /// let reader = bgzf::Reader::with_gzi(&data[..], gzi);
    /// assert_eq!(reader.uncompressed_position(), Some(0));
    /// ```
    pub fn uncompressed_position(&self) -> Option<u64> {
        self.uncompressed_position
    }

    /// Returns the current uncompressed position of the stream.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bgzf as bgzf;
    /// let data = [];
    /// let gzi = vec![];
    /// let reader = bgzf::Reader::with_gzi(&data[..], gzi);
    /// assert_eq!(reader.get_gzi(), Some(&vec![]));
    /// ```
    pub fn get_gzi(&self) -> Option<&gzi::Index> {
        self.gzi.as_ref()
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
    /// reader.seek_virtual_position(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek_virtual_position(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        let (cpos, upos) = pos.into();

        self.inner.get_mut().seek(SeekFrom::Start(cpos))?;
        self.position = cpos;

        self.read_block()?;

        self.block.data_mut().set_position(usize::from(upos));

        if let Some(gzi) = self.get_gzi() {
            let gzi = [&[(0, 0)], &gzi[..]].concat();
            if let Ok(i) = gzi.binary_search_by(|p| p.0.cmp(&cpos)) {
                self.uncompressed_position = Some(gzi[i].1 + upos as u64);
            } else {
                return Err(io::Error::new(io::ErrorKind::Other, "fail to find cpos in gzi"));
            }
        }

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
        let prev = self.block.data().position();
        self.block.data_mut().consume(amt);
        if let Some(pos) = self.uncompressed_position() {
            let offset = self.block.data().position() - prev;
            let pos = pos
                .checked_add(offset as u64)
                .expect("overflow when adding offset to uncompressed position");
            self.uncompressed_position = Some(pos);
        }
    }

    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.block.data().has_remaining() {
            self.read_block()?;
        }

        Ok(self.block.data().as_ref())
    }
}

impl<R> Seek for Reader<R>
where
    R: Read + Seek,
{
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        use io::{Error, ErrorKind};

        if let Reader {
            gzi: Some(gzi),
            uncompressed_position: Some(cur_upos),
            ..
        } = self {
            let gzi = [&[(0, 0)], &gzi[..]].concat();
            let cur_upos = *cur_upos;
            let cur_vpos = self.virtual_position();
            let mut is_seek_end = false;
            let new_pos = match pos {
                SeekFrom::Current(n) => {
                    if n < 0 {
                        cur_upos.checked_sub((-n) as u64)
                    } else {
                        cur_upos.checked_add(n as u64)
                    }
                    .expect("overflow when adding offset to uncompressed position")
                }
                SeekFrom::Start(n) => n,
                SeekFrom::End(n) => {
                    is_seek_end = true;
                    let end_cpos = gzi.last().map(|p| p.0).unwrap();
                    if let Ok(end_pos) = VirtualPosition::try_from((end_cpos, 0)) {
                        self.seek_virtual_position(end_pos)?;
                        let end_upos = self.uncompressed_position().unwrap()
                            + self.block.data().len() as u64;
                        if n < 0 {
                            end_upos.checked_sub((-n) as u64)
                        } else {
                            end_upos.checked_add(n as u64)
                        }
                        .expect("overflow when adding offset to the end of uncompressed position")
                    } else {
                        return Err(Error::new(ErrorKind::Other, "fail to convert tuple to VirtualPosition"));
                    }
                }
            };
            let i = gzi.partition_point(|p| p.1 <= new_pos);
            let (cpos, gzi_upos) = gzi[i - 1];
            let upos = new_pos - gzi_upos;
            if upos > u16::MAX as u64 {
                if is_seek_end {
                    self.seek_virtual_position(cur_vpos)?;
                }
                return Err(Error::new(ErrorKind::Other, ""))
            } else {
                if let Ok(virtual_position) = VirtualPosition::try_from((cpos, upos as u16)) {
                    self.seek_virtual_position(virtual_position)?;
                } else {
                    if is_seek_end {
                        self.seek_virtual_position(cur_vpos)?;
                    }
                    return Err(Error::new(ErrorKind::Other, "fail to convert tuple to VirtualPosition")); 
                }
            }
            Ok(new_pos)
        } else {
            return Err(Error::new(ErrorKind::NotFound, "index not found"));
        }
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
    fn test_seek_virtual_position() -> Result<(), Box<dyn std::error::Error>> {
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

        reader.seek_virtual_position(VirtualPosition::try_from((0, 3))?)?;

        buf.clear();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"dles");
        assert_eq!(reader.virtual_position(), eof);

        Ok(())
    }

    #[test]
    fn test_seek() -> Result<(), Box<dyn std::error::Error>> {
        #[rustfmt::skip]
        let data = [
            // block 0
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x29, 0x00, 0xb3, 0xcb, 0xcb, 0xcf, 0x4f, 0xc9, 0x49, 0x2d, 0xe6, 0x72,
            0x0c, 0x71, 0x76, 0xe7, 0x02, 0x00, 0x66, 0x54, 0x8f, 0x56, 0x0e, 0x00, 0x00, 0x00,
            // EOF block
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        let gzi = [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];

        let mut gzi_reader = gzi::Reader::new(&gzi[..]);
        let gzi = gzi_reader.read_index()?;

        let eof = 14;

        let mut reader =
            Reader::with_gzi(Cursor::new(&data), gzi);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(reader.uncompressed_position().unwrap(), eof);

        reader.seek(SeekFrom::Start(1))?;

        let mut line = String::new();
        reader.read_line(&mut line)?;
        assert_eq!(line, "noodles\n");

        reader.seek(SeekFrom::End(-13))?;
        line.clear();
        reader.read_line(&mut line)?;
        assert_eq!(line, "noodles\n");

        reader.seek(SeekFrom::Current(2))?;
        buf.clear();
        reader.read_to_end(&mut buf)?;
        assert_eq!(buf, b"CG\n");

        assert_eq!(reader.uncompressed_position().unwrap(), eof);

        Ok(())
    }
}
