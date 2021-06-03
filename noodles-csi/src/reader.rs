use std::{
    convert::TryFrom,
    io::{self, Read},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf::{self as bgzf, index::Chunk};

use super::{
    index::{
        reference_sequence::{Bin, Metadata},
        ReferenceSequence,
    },
    Index, MAGIC_NUMBER,
};

/// A CSI reader.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a CSI reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(inner: R) -> Self {
        Self {
            inner: bgzf::Reader::new(inner),
        }
    }

    /// Reads a CSI index.
    ///
    /// The position of the stream is expected to be at the beginning.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_csi as csi;
    /// let mut reader = File::open("sample.bcf.csi").map(csi::Reader::new)?;
    /// let index = reader.read_index();
    /// # Ok::<(), io::Error>(())
    /// ```

    pub fn read_index(&mut self) -> io::Result<Index> {
        read_magic(&mut self.inner)?;

        let min_shift = self.inner.read_i32::<LittleEndian>()?;
        let depth = self.inner.read_i32::<LittleEndian>()?;
        let aux = read_aux(&mut self.inner)?;
        let reference_sequences = read_reference_sequences(&mut self.inner, depth)?;
        let n_no_coor = self.inner.read_u64::<LittleEndian>().ok();

        let mut builder = Index::builder()
            .set_min_shift(min_shift)
            .set_depth(depth)
            .set_aux(aux)
            .set_reference_sequences(reference_sequences);

        if let Some(n_no_coor) = n_no_coor {
            builder = builder.set_n_no_coor(n_no_coor);
        }

        Ok(builder.build())
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;

    if magic == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid CSI file format",
        ))
    }
}

fn read_aux<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let l_aux = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut aux = vec![0; l_aux];
    reader.read_exact(&mut aux)?;

    Ok(aux)
}

fn read_reference_sequences<R>(reader: &mut R, depth: i32) -> io::Result<Vec<ReferenceSequence>>
where
    R: Read,
{
    let n_ref = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = Vec::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (bins, metadata) = read_bins(reader, depth)?;
        let reference_sequence = ReferenceSequence::new(bins, metadata);
        reference_sequences.push(reference_sequence);
    }

    Ok(reference_sequences)
}

fn read_bins<R>(reader: &mut R, depth: i32) -> io::Result<(Vec<Bin>, Option<Metadata>)>
where
    R: Read,
{
    let n_bin = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut bins = Vec::with_capacity(n_bin);

    let metadata_bin_id = Bin::metadata_id(depth);
    let mut metadata = None;

    for _ in 0..n_bin {
        let id = reader.read_u32::<LittleEndian>()?;
        let loffset = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;
        let chunks = read_chunks(reader)?;

        let bin = Bin::new(id, loffset, chunks);

        if id == metadata_bin_id {
            metadata = Metadata::try_from(&bin)
                .map(Some)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        } else {
            bins.push(bin);
        }
    }

    Ok((bins, metadata))
}

fn read_chunks<R>(reader: &mut R) -> io::Result<Vec<Chunk>>
where
    R: Read,
{
    let n_chunk = reader.read_i32::<LittleEndian>().and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut chunks = Vec::with_capacity(n_chunk);

    for _ in 0..n_chunk {
        let chunk_beg = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        let chunk_end = reader
            .read_u64::<LittleEndian>()
            .map(bgzf::VirtualPosition::from)?;

        chunks.push(Chunk::new(chunk_beg, chunk_end));
    }

    Ok(chunks)
}
