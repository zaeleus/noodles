pub use self::{query::Query, records::Records, references::References};

mod query;
mod records;
mod references;

use std::{
    ffi::CStr,
    io::{self, Read, Seek},
    ops::DerefMut,
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles::formats::bai;
use noodles_bgzf as bgzf;

use super::{Record, Reference, MAGIC_NUMBER};

pub type Header = String;
pub type Meta = (Header, Vec<Reference>);

const BLOCK_SIZE_LEN: usize = 4;

pub struct Reader<R: Read> {
    inner: bgzf::Reader<R>,
    block: bgzf::Block,
}

impl<R: Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
            block: bgzf::Block::default(),
        }
    }

    #[deprecated(note = "Use reader.header() instead.")]
    pub fn meta(&mut self) -> io::Result<Meta> {
        self.header()
    }

    pub fn header(&mut self) -> io::Result<Meta> {
        self.inner.read_block(&mut self.block)?;

        let mut reader = self.block.deref_mut();

        let magic = read_magic(&mut reader)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM header",
            ));
        }

        let header = read_header(&mut reader)?;
        let references = read_references(&mut reader)?;

        Ok((header, references))
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let mut buf = [0; BLOCK_SIZE_LEN];

        if let Err(e) = self.read_exact(&mut buf) {
            match e.kind() {
                io::ErrorKind::UnexpectedEof => return Ok(0),
                _ => return Err(e),
            }
        }

        let block_size = u32::from_le_bytes(buf) as usize;

        record.resize(block_size);
        let buf = record.deref_mut();

        self.read_exact(buf)?;

        Ok(block_size)
    }

    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
    }

    pub fn virtual_position(&self) -> u64 {
        self.block.virtual_position()
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        let mut bytes_read = 0;

        loop {
            match self.block.read(&mut buf[bytes_read..]) {
                Ok(0) => match self.inner.read_block(&mut self.block) {
                    Ok(0) => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
                    Ok(_) => {}
                    Err(e) => return Err(e),
                },
                Ok(len) => {
                    bytes_read += len;

                    if bytes_read == buf.len() {
                        break;
                    }
                }
                Err(e) => return Err(e),
            }
        }

        Ok(())
    }
}

impl<R: Read + Seek> Reader<R> {
    pub fn seek(&mut self, pos: u64) -> io::Result<u64> {
        self.inner.seek(pos, &mut self.block)
    }

    pub fn query(&mut self, index_ref: &bai::Reference, start: i32, end: i32) -> Query<R> {
        let query_bins = bai::query(&index_ref.bins, start as u64, end as u64);

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .cloned()
            .collect();

        let merged_chunks = bai::merge_chunks(&chunks);

        Query::new(self, merged_chunks)
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
