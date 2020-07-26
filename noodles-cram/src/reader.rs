pub mod block;
mod compression_header;
mod container;
mod encoding;
pub mod record;
mod slice;

use std::{
    io::{self, Read},
    str,
};

use byteorder::{LittleEndian, ReadBytesExt};

use crate::{container::Slice, Container};

use super::MAGIC_NUMBER;

use self::{block::read_block, compression_header::read_compression_header};

pub struct Reader<R>
where
    R: Read,
{
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub fn new(reader: R) -> Self {
        Self { inner: reader }
    }

    pub fn read_file_definition(&mut self) -> io::Result<[u8; 20]> {
        let magic = read_magic(&mut self.inner)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid CRAM header",
            ));
        }

        read_format(&mut self.inner)?;

        read_file_id(&mut self.inner)
    }

    pub fn read_file_header(&mut self) -> io::Result<String> {
        let container_header = container::read_header(&mut self.inner)?;

        let mut buf = vec![0; container_header.len() as usize];
        self.inner.read_exact(&mut buf)?;

        let mut reader = &buf[..];
        let block = read_block(&mut reader)?;

        let mut data = &block.decompressed_data()[..];
        let _header_len = data.read_i32::<LittleEndian>()?;

        str::from_utf8(data)
            .map(|s| s.into())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub fn read_container(&mut self, container: &mut Container) -> io::Result<()> {
        let header = container::read_header(&mut self.inner)?;
        let compression_header = read_compression_header(&mut self.inner)?;

        container.slices_mut().clear();

        for _ in 0..header.landmarks().len() {
            let slice_header = slice::read_header(&mut self.inner)?;
            let mut slice = Slice::new(slice_header);
            slice::read_blocks(&mut self.inner, &mut slice)?;
            container.add_slice(slice);
        }

        *container.header_mut() = header;
        *container.compression_header_mut() = compression_header;

        Ok(())
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut buf = [0; 4];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

fn read_format<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: Read,
{
    let mut buf = [0; 2];
    reader.read_exact(&mut buf)?;
    Ok((buf[0], buf[1]))
}

fn read_file_id<R>(reader: &mut R) -> io::Result<[u8; 20]>
where
    R: Read,
{
    let mut buf = [0; 20];
    reader.read_exact(&mut buf)?;
    Ok(buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_file_definition() -> Vec<u8> {
        let mut data = MAGIC_NUMBER.to_vec();

        let format = [0x03, 0x00];
        data.extend_from_slice(&format);

        let file_id = [
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f,
        ];
        data.extend_from_slice(&file_id);

        data
    }

    #[test]
    fn test_read_file_definition() -> io::Result<()> {
        let data = build_file_definition();
        let mut reader = Reader::new(&data[..]);
        let file_id = reader.read_file_definition()?;

        let expected = [
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f,
        ];

        assert_eq!(file_id, expected);

        Ok(())
    }

    #[test]
    fn test_read_file_definition_with_invalid_magic_number() {
        let data = b"BAM\x01";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_definition().is_err());
    }
}
