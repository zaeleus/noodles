pub mod block;
mod container;
mod slice;

use std::io::{self, Read};

use crate::Container;

static MAGIC_NUMBER: &[u8] = b"CRAM";

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

    pub fn read_container(&mut self, container: &mut Container) -> io::Result<()> {
        let header = container::read_header(&mut self.inner)?;

        let blocks = container.blocks_mut();
        blocks.resize(header.len() as usize, Default::default());
        self.inner.read_exact(blocks)?;

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
    fn test_read_file_definition() {
        let data = build_file_definition();
        let mut reader = Reader::new(&data[..]);
        let file_id = reader.read_file_definition().unwrap();

        let expected = [
            0x00, 0x68, 0xac, 0xf3, 0x06, 0x4d, 0xaa, 0x1e, 0x29, 0xa4, 0xa0, 0x8c, 0x56, 0xee,
            0x91, 0x9b, 0x91, 0x04, 0x21, 0x1f,
        ];

        assert_eq!(file_id, expected);
    }

    #[test]
    fn test_read_file_definition_with_invalid_magic_number() {
        let data = b"BAM\x01";
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_definition().is_err());
    }
}
