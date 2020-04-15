pub use self::{block::Block, reader::Reader, virtual_position::VirtualPosition, writer::Writer};

mod block;
mod gz;
mod reader;
mod virtual_position;
mod writer;

const GZIP_XLEN_SIZE: usize = 2;
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
        writer.flush()?;

        let data = writer.get_ref();
        let mut reader = Reader::new(&data[..]);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf)?;

        assert_eq!(buf, b"noodles-bgzf");

        Ok(())
    }
}
