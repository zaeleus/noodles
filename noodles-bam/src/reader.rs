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
use noodles::Region;
use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_sam::header::ReferenceSequence;

use super::{bai, Record, Reference, MAGIC_NUMBER};

const BLOCK_SIZE_LEN: usize = 4;

pub struct Reader<R: Read> {
    inner: bgzf::Reader<R>,
}

impl<R: Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    pub fn read_header(&mut self) -> io::Result<String> {
        let magic = read_magic(&mut self.inner)?;

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM header",
            ));
        }

        let header = read_header(&mut self.inner)?;

        // I'm actually not sure what the BAM references are used for since the reference sequence
        // dictionary in the SAM header should have or has the same information plus, commonly,
        // more (e.g., the M5 tag).
        //
        // In this case, the BAM references are discarded in favor of the dictionary in the SAM
        // header.
        read_references(&mut self.inner)?;

        Ok(header)
    }

    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let mut buf = [0; BLOCK_SIZE_LEN];

        if let Err(e) = self.inner.read_exact(&mut buf) {
            match e.kind() {
                io::ErrorKind::UnexpectedEof => return Ok(0),
                _ => return Err(e),
            }
        }

        let block_size = u32::from_le_bytes(buf) as usize;

        record.resize(block_size);
        let buf = record.deref_mut();

        self.inner.read_exact(buf)?;

        Ok(block_size)
    }

    pub fn records(&mut self) -> Records<R> {
        Records::new(self)
    }

    pub fn virtual_position(&self) -> VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R: Read + Seek> Reader<R> {
    pub fn seek(&mut self, pos: VirtualPosition) -> io::Result<VirtualPosition> {
        self.inner.seek(pos)
    }

    pub fn query(
        &mut self,
        reference_sequences: &[ReferenceSequence],
        index: &bai::Index,
        region: &Region,
    ) -> io::Result<Query<R>> {
        let (i, _, start, end) = region.resolve(reference_sequences).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("could not resolve region: {:?}", region),
            )
        })?;

        let index_reference = index.references().get(i).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "could not find reference in index: {} >= {}",
                    i,
                    reference_sequences.len()
                ),
            )
        })?;

        let query_bins = bai::query(index_reference.bins(), start, end);

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .cloned()
            .collect();

        let merged_chunks = bai::merge_chunks(&chunks);

        Ok(Query::new(self, merged_chunks, i, start, end))
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

fn read_header<R>(reader: &mut R) -> io::Result<String>
where
    R: Read,
{
    let l_text = reader.read_i32::<LittleEndian>()?;

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
    let n_ref = reader.read_i32::<LittleEndian>()?;
    References::new(reader, n_ref as usize).collect()
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
