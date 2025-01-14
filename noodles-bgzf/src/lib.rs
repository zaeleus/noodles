#![warn(missing_docs)]

//! **noodles-bgzf** handles the reading and writing of the blocked gzip format (BGZF).
//!
//! While the gzip format is typically a single stream, a BGZF is the concatenation of many gzip
//! streams. Each stream is called a block, with its uncompressed data size being constrained to
//! less than 64 KiB. This multistream gzip allows random access using [`virtual positions`].
//!
//! noodles-bgzf abstracts away the concept of blocks, implementing [`std::io::Read`] for the
//! reader and [`std::io::Write`] for the writer.
//!
//! [`virtual positions`]: VirtualPosition
//!
//! # Examples
//!
//! ## Read an entire BGZF file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, Read}};
//! use noodles_bgzf as bgzf;
//! let mut reader = File::open("data.gz").map(bgzf::Reader::new)?;
//! let mut data = Vec::new();
//! reader.read_to_end(&mut data)?;
//! # Ok::<(), io::Error>(())
//! ```
//!
//! ## Write a BGZF file
//!
//! ```no_run
//! # use std::{fs::File, io::{self, Write}};
//! use noodles_bgzf as bgzf;
//! let mut writer = File::create("data.gz").map(bgzf::Writer::new)?;
//! writer.write_all(b"noodles-bgzf")?;
//! # Ok::<(), io::Error>(())
//! ```

#[cfg(feature = "async")]
pub mod r#async;

mod block;
pub(crate) mod deflate;
mod gz;
pub mod gzi;
pub mod indexed_reader;
pub mod io;
mod multithreaded_reader;
pub mod multithreaded_writer;
pub mod reader;
pub mod virtual_position;
pub mod writer;

pub use self::{
    indexed_reader::IndexedReader, multithreaded_reader::MultithreadedReader,
    multithreaded_writer::MultithreadedWriter, reader::Reader, virtual_position::VirtualPosition,
    writer::Writer,
};

#[cfg(feature = "async")]
#[deprecated(since = "0.35.0", note = "Use `bgzf::r#async::Reader` instead.")]
pub use self::r#async::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.35.0", note = "Use `bgzf::r#async::Writer` instead.")]
pub use self::r#async::Writer as AsyncWriter;

use self::block::Block;

// XLEN (2)
const GZIP_XLEN_SIZE: usize = 2;

// SI1 (1) + SI2 (1) + SLEN (2) + BSIZE (2)
const BGZF_XLEN: usize = 6;

// § 4.1 The BGZF compression format (2021-06-03): "Thus while `ISIZE` is stored as a `uint32_t` as
// per the gzip format, in BGZF it is limited to the range [0, 65536]."
const BGZF_MAX_ISIZE: usize = 1 << 16;

pub(crate) const BGZF_HEADER_SIZE: usize = gz::HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;

#[cfg(test)]
mod tests {
    use std::io::{self, BufRead, Cursor, Read, Write};

    use super::*;

    #[test]
    fn test_self() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        writer.write_all(b"noodles")?;
        writer.flush()?;
        writer.write_all(b"-")?;
        writer.flush()?;
        writer.write_all(b"bgzf")?;

        let data = writer.finish()?;
        let mut reader = Reader::new(&data[..]);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"noodles-bgzf");

        Ok(())
    }

    #[test]
    fn test_self_buffered() -> io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        writer.write_all(b"noodles\n-\nbgzf\nbuffered")?;

        let data = writer.finish()?;
        let mut reader = Reader::new(&data[..]);

        let mut lines = Vec::new();
        let mut virtual_positions = Vec::new();

        loop {
            virtual_positions.push(reader.virtual_position());

            let mut line = String::new();
            match reader.read_line(&mut line) {
                Ok(0) => {
                    virtual_positions.pop();
                    break;
                }
                Err(e) => return Err(e),
                _ => (),
            }

            lines.push(line);
        }

        let expected_lines = vec!["noodles\n", "-\n", "bgzf\n", "buffered"];
        assert_eq!(lines, expected_lines);

        let expected_upos = [0, 8, 10, 15];
        let expected_virtual_positions: Vec<VirtualPosition> = expected_upos
            .iter()
            .map(|x| VirtualPosition::try_from((0, *x)).unwrap())
            .collect();
        assert_eq!(virtual_positions, expected_virtual_positions);

        Ok(())
    }

    #[test]
    fn test_self_multithreaded() -> io::Result<()> {
        let mut writer = MultithreadedWriter::new(Vec::new());

        writer.write_all(b"noodles")?;
        writer.flush()?;
        writer.write_all(b"-")?;
        writer.flush()?;
        writer.write_all(b"bgzf")?;

        let data = writer.finish().map(Cursor::new)?;
        let mut reader = MultithreadedReader::new(data);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"noodles-bgzf");

        Ok(())
    }
}
