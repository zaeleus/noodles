#![deny(missing_docs)]

//! **noodles-bgzf** handles the reading and writing of BGZF (blocked gzip file).
//!
//! While the gzip format is typically a single stream, a BGZF is the concatentation of many gzip
//! streams. Each stream is called a block, with its uncompressed data size being constrained to
//! less than 64 KiB. This multistream gzip allows random access using [`virtual positions`].
//!
//! noodles-bgzf abstracts away the concept of blocks, implementing [`std::io::Read`] for the
//! reader and [`std::io::Write`] for the writer.
//!
//! [`virtual positions`]: struct.VirtualPosition.html
//! [`std::io::Read`]: https://doc.rust-lang.org/std/io/trait.Read.html
//! [`std::io::Write`]: https://doc.rust-lang.org/std/io/trait.Write.html
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
//! reader.read_to_end(&mut data);
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

pub use self::{block::Block, reader::Reader, virtual_position::VirtualPosition, writer::Writer};

mod block;
mod gz;
mod reader;
pub mod virtual_position;
mod writer;

// XLEN (2)
const GZIP_XLEN_SIZE: usize = 2;

// SI1 (1) + SI2 (1) + SLEN (2) + BSIZE (2)
const BGZF_XLEN: usize = 6;

pub(crate) const BGZF_HEADER_SIZE: usize = gz::HEADER_SIZE + GZIP_XLEN_SIZE + BGZF_XLEN;

#[cfg(test)]
mod tests {
    use std::io::{self, Read, Write};

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
}
