use std::{
    ffi::CStr,
    io::{self, Read, Seek},
    ops::DerefMut,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles::formats::bgzf;

use super::{Record, Reference, MAGIC_NUMBER};

pub type Header = String;
pub type Meta = (Header, Vec<Reference>, bgzf::Block);

pub struct Reader2<R: Read + Seek> {
    inner: bgzf::Reader<R>,
}

impl<R: Read + Seek> Reader2<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    pub fn meta(&mut self) -> io::Result<Meta> {
        let mut block = bgzf::Block::new();
        self.inner.read_block(&mut block)?;

        let mut reader = block.deref_mut();

        let magic = read_magic(&mut reader)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM header",
            ));
        }

        let header = read_header(&mut reader)?;
        let references = read_references(&mut reader)?;

        Ok((header, references, block))
    }

    pub fn read_record(
        &mut self,
        block: &mut bgzf::Block,
        record: &mut Record,
    ) -> io::Result<usize> {
        if block.is_eof() {
            self.inner.read_block(block)?;
        }

        let block_size = match block.read_block_size() {
            Ok(bs) => bs as usize,
            Err(_) => return Ok(0),
        };

        record.resize(block_size);
        let mut buf = record.deref_mut();

        while !buf.is_empty() {
            match block.read_record(buf) {
                Ok(0) => break,
                Ok(n) => {
                    if block.is_eof() {
                        match self.inner.read_block(block) {
                            Ok(0) => return Ok(0),
                            Ok(_) => {}
                            Err(e) => return Err(e),
                        }
                    }

                    let tmp = buf;
                    buf = &mut tmp[n..];
                }
                Err(e) => return Err(e),
            }
        }

        Ok(block_size)
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<[u8; 4]>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    Ok(magic)
}

fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: Read,
{
    let l_text = read_l_text(reader)?;

    let mut c_text = vec![0; l_text as usize];
    reader.read_exact(&mut c_text)?;

    // Headers are not necessarily NUL-terminated.
    bytes_with_nul_to_string(&c_text).or_else(|_| {
        String::from_utf8(c_text).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_references<R>(reader: &mut R) -> io::Result<Vec<Reference>>
where
    R: Read,
{
    let n_ref = read_n_ref(reader)?;

    let references = References::new(reader, n_ref as usize)
        .map(|r| r.unwrap())
        .collect();

    Ok(references)
}

fn read_l_text<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

fn read_n_ref<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

struct References<'a, R: Read> {
    reader: &'a mut R,
    i: usize,
    len: usize,
}

impl<'a, R: Read> References<'a, R> {
    fn new(reader: &'a mut R, len: usize) -> References<R> {
        References { reader, i: 0, len }
    }
}

impl<'a, R: 'a + Read> Iterator for References<'a, R> {
    type Item = io::Result<Reference>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i >= self.len {
            return None;
        }

        let result = read_reference(self.reader);

        self.i += 1;

        Some(result)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.len, Some(self.len))
    }
}

fn read_reference<R>(reader: &mut R) -> io::Result<Reference>
where
    R: Read,
{
    let l_name = read_l_name(reader)?;
    let name = read_name(reader, l_name as usize)?;
    let l_ref = read_l_ref(reader)?;
    Ok(Reference::new(name, l_ref))
}

fn read_l_name<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

fn read_name<R>(reader: &mut R, l_name: usize) -> io::Result<String>
where
    R: Read,
{
    let mut buf = vec![0; l_name];
    reader.read_exact(&mut buf)?;
    bytes_with_nul_to_string(&buf)
}

fn read_l_ref<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    reader.read_i32::<LittleEndian>()
}

fn bytes_with_nul_to_string(buf: &[u8]) -> io::Result<String> {
    CStr::from_bytes_with_nul(buf)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|c_str| {
            c_str
                .to_str()
                .map(|s| s.to_string())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}
