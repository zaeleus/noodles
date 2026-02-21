//! BGZF I/O.

mod block;
mod buf_read;
pub mod indexed_reader;
mod multithreaded_reader;
pub mod multithreaded_writer;
mod read;
pub mod reader;
mod seek;
pub mod writer;

pub(crate) use self::block::Block;
pub use self::{
    buf_read::BufRead,
    indexed_reader::IndexedReader,
    multithreaded_reader::MultithreadedReader,
    multithreaded_writer::{BlockInfo, BlockInfoRx, MultithreadedWriter},
    read::Read,
    reader::Reader,
    seek::Seek,
    writer::Writer,
};

#[cfg(test)]
mod tests {
    use std::io::{self, BufRead, Cursor, Read, Write};

    use super::*;
    use crate::VirtualPosition;

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

    #[test]
    fn test_multithreaded_position_tracking() -> io::Result<()> {
        use std::num::NonZero;

        let mut writer = multithreaded_writer::Builder::default()
            .set_worker_count(NonZero::new(2).unwrap())
            .build_from_writer(Vec::new());

        // Get the receiver before writing
        let rx = writer.block_info_receiver().unwrap().clone();

        // Initial state
        assert_eq!(writer.position(), 0);
        assert_eq!(writer.current_block_number(), 0);
        assert_eq!(writer.blocks_written(), 0);
        assert_eq!(writer.buffer_offset(), 0);

        // Write first block
        writer.write_all(b"hello")?;
        assert_eq!(writer.buffer_offset(), 5);
        writer.flush()?;
        assert_eq!(writer.current_block_number(), 1);
        assert_eq!(writer.buffer_offset(), 0);

        // Write second block
        writer.write_all(b"world")?;
        writer.flush()?;
        assert_eq!(writer.current_block_number(), 2);

        // Finish to ensure all blocks are written
        let _data = writer.finish()?;

        // Collect block info
        let mut block_infos: Vec<_> = rx.try_iter().collect();
        block_infos.sort_by_key(|b| b.block_number);

        // Should have 2 blocks
        assert_eq!(block_infos.len(), 2);

        // First block
        assert_eq!(block_infos[0].block_number, 0);
        assert_eq!(block_infos[0].compressed_start, 0);
        assert_eq!(block_infos[0].uncompressed_size, 5);
        assert!(block_infos[0].compressed_size > 0);

        // Second block starts after first
        assert_eq!(block_infos[1].block_number, 1);
        assert_eq!(
            block_infos[1].compressed_start,
            block_infos[0].compressed_size as u64
        );
        assert_eq!(block_infos[1].uncompressed_size, 5);

        Ok(())
    }
}
