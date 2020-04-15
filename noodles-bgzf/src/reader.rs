use std::io::{self, Read, Seek, SeekFrom};

use byteorder::{ByteOrder, LittleEndian};
use flate2::read::DeflateDecoder;

use super::{gz, Block, VirtualPosition, BGZF_HEADER_SIZE};

pub struct Reader<R: Read> {
    inner: R,
    position: u64,
    cdata: Vec<u8>,
    block: Block,
}

impl<R: Read> Reader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            position: 0,
            cdata: Vec::new(),
            block: Block::default(),
        }
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn read_block(&mut self, block: &mut Block) -> io::Result<usize> {
        match read_block(&mut self.inner, &mut self.cdata, block) {
            Ok(0) => Ok(0),
            Ok(bs) => {
                block.set_c_offset(self.position);
                self.position += bs as u64;
                Ok(bs)
            }
            Err(e) => Err(e),
        }
    }
}

impl<R: Read + Seek> Reader<R> {
    pub fn seek(&mut self, pos: VirtualPosition, block: &mut Block) -> io::Result<VirtualPosition> {
        let compressed_offset = pos.compressed();
        let uncompressed_offset = pos.uncompressed();

        self.inner.seek(SeekFrom::Start(compressed_offset))?;
        self.position = compressed_offset;

        self.read_block(block)?;
        block.seek(SeekFrom::Start(uncompressed_offset))?;

        Ok(pos)
    }
}

impl<R> Read for Reader<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self.block.read(buf) {
            Ok(0) => match read_block(&mut self.inner, &mut self.cdata, &mut self.block) {
                Ok(0) => Ok(0),
                Ok(bs) => {
                    self.block.set_c_offset(self.position);
                    self.position += bs as u64;
                    Err(io::Error::from(io::ErrorKind::Interrupted))
                }
                Err(e) => Err(e),
            },
            Ok(n) => Ok(n),
            Err(e) => Err(e),
        }
    }
}

fn read_block_size<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut header = [0; BGZF_HEADER_SIZE];

    if reader.read_exact(&mut header).is_err() {
        return Ok(0);
    }

    let bsize = &header[16..18];

    // Add 1 because BSIZE is "total Block SIZE minus 1".
    Ok(LittleEndian::read_u16(bsize) + 1)
}

fn read_trailer<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut trailer = [0; gz::TRAILER_SIZE];
    reader.read_exact(&mut trailer)
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
    let block_size = match read_block_size(reader).map(usize::from) {
        Ok(0) => return Ok(0),
        Ok(bs) => bs,
        Err(e) => return Err(e),
    };

    let cdata_len = block_size - BGZF_HEADER_SIZE - gz::TRAILER_SIZE;
    cdata.resize(cdata_len, Default::default());
    reader.read_exact(cdata)?;

    read_trailer(reader)?;

    let block_buf = block.get_mut();
    block_buf.clear();

    inflate_data(&cdata[..], block_buf)?;

    block.set_position(0);

    Ok(block_size)
}

#[cfg(test)]
mod tests {
    use std::fs::File;

    use super::*;

    #[test]
    fn test_read_block() -> io::Result<()> {
        let file = File::open("tests/fixtures/sample.gz")?;
        let mut reader = Reader::new(file);

        let mut block = Block::default();

        reader.read_block(&mut block)?;
        assert_eq!(&block.get_ref()[..], &b"noodles"[..]);

        reader.read_block(&mut block)?;
        assert_eq!(&block.get_ref()[..], &b"bgzf"[..]);

        reader.read_block(&mut block)?;
        assert!(block.get_ref().is_empty());

        assert_eq!(reader.read_block(&mut block)?, 0);

        Ok(())
    }
}
