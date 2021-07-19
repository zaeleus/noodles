mod query;
pub mod record;
mod records;
pub mod string_map;
pub mod value;

pub use self::{query::Query, records::Records};

use std::{
    convert::TryFrom,
    ffi::CStr,
    io::{self, Read, Seek},
};

use byteorder::{LittleEndian, ReadBytesExt};
use noodles_bgzf as bgzf;
use noodles_core::{region::Interval, Region};
use noodles_csi::{BinningIndex, BinningIndexReferenceSequence};
use noodles_vcf::header::Contigs;

use super::{Record, MAGIC_NUMBER};

/// A BCF reader.
///
/// The BCF format is comprised of two parts: 1) a VCF header and 2) a list of records.
pub struct Reader<R> {
    inner: bgzf::Reader<R>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BCF reader.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn new(reader: R) -> Self {
        Self {
            inner: bgzf::Reader::new(reader),
        }
    }

    /// Reads the BCF file format.
    ///
    /// The BCF magic number is also checked.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the major and minor format versions as a tuple.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// let (major, minor) = reader.read_file_format()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_file_format(&mut self) -> io::Result<(u8, u8)> {
        let mut buf = [0; 5];

        self.inner.read_exact(&mut buf)?;

        let magic = &buf[..3];

        if magic != MAGIC_NUMBER {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid file format",
            ));
        }

        let major = buf[3];
        let minor = buf[4];

        Ok((major, minor))
    }

    /// Reads the raw VCF header.
    ///
    /// The position of the stream is expected to be directly after the file format.
    ///
    /// This returns the raw VCF header as a [`std::string::String`]. It can subsequently be parsed
    /// as a [`noodles_vcf::Header`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// reader.read_file_format()?;
    /// let header = reader.read_header()?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        let l_text = self.inner.read_u32::<LittleEndian>().and_then(|n| {
            usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })?;

        let mut buf = vec![0; l_text];
        self.inner.read_exact(&mut buf)?;

        CStr::from_bytes_with_nul(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|c_header| {
                c_header
                    .to_str()
                    .map(|s| s.into())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }

    /// Reads a single record.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// It is more ergnomic to read records using an iterator (see [`Self::records`]), but using
    /// this method directly allows the reuse of a single [`Record`] buffer.
    ///
    /// If successful, the record size is returned. If a record size of 0 is returned, the stream
    /// reached EOF.
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        let l_shared = match self.inner.read_u32::<LittleEndian>() {
            Ok(len) => len,
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
            Err(e) => return Err(e),
        };

        let l_indiv = self.inner.read_u32::<LittleEndian>()?;

        let record_len = l_shared
            .checked_add(l_indiv)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "invalid record length: l_shared = {}, l_indiv = {}",
                        l_shared, l_indiv
                    ),
                )
            })
            .and_then(|len| {
                usize::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

        record.resize(record_len);
        self.inner.read_exact(record)?;

        Ok(record_len)
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be directly after the header or at the start of another record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// reader.read_file_format()?;
    /// reader.read_header()?;
    ///
    /// for result in reader.records() {
    ///     let record = result?;
    ///     println!("{:?}", &record[..]);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }

    /// Returns the current virtual position of the underlying BGZF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    ///
    /// let data = Vec::new();
    /// let reader = bcf::Reader::new(&data[..]);
    /// let virtual_position = reader.virtual_position();
    ///
    /// assert_eq!(virtual_position.compressed(), 0);
    /// assert_eq!(virtual_position.uncompressed(), 0);
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn virtual_position(&self) -> bgzf::VirtualPosition {
        self.inner.virtual_position()
    }
}

impl<R> Reader<R>
where
    R: Read + Seek,
{
    /// Seeks the underlying BGZF reader to the given virtual position.
    ///
    /// Virtual positions typically come from an associated BCF index file.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_bcf as bcf;
    /// use noodles_bgzf as bgzf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    ///
    /// let virtual_position = bgzf::VirtualPosition::from(102334155);
    /// reader.seek(virtual_position)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn seek(&mut self, pos: bgzf::VirtualPosition) -> io::Result<bgzf::VirtualPosition> {
        self.inner.seek(pos)
    }

    /// Returns an iterator over records that intersects the given region.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    /// use noodles_core::Region;
    /// use noodles_csi as csi;
    /// use noodles_vcf as vcf;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    ///
    /// let header: vcf::Header = reader.read_header()?.parse()?;
    /// let contigs = header.contigs();
    ///
    /// let index = csi::read("sample.bcf.csi")?;
    /// let region = Region::mapped("sq0", 8..=13);
    /// let query = reader.query(&contigs, &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query<I, IRS>(
        &mut self,
        contigs: &Contigs,
        index: &I,
        region: &Region,
    ) -> io::Result<Query<'_, R>>
    where
        I: BinningIndex<IRS>,
        IRS: BinningIndexReferenceSequence,
    {
        let (reference_sequence_id, interval) = resolve_region(contigs, region)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        Ok(Query::new(self, chunks, reference_sequence_id, interval))
    }
}

fn resolve_region(contigs: &Contigs, region: &Region) -> io::Result<(usize, Interval)> {
    if let Some(r) = region.as_mapped() {
        let i = contigs.get_index_of(r.name()).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("region does not exist in contigs: {:?}", region),
            )
        })?;

        Ok((i, r.interval()))
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "region is not mapped",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn compress(data: &[u8]) -> io::Result<Vec<u8>> {
        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(data)?;
        writer.finish()
    }

    #[test]
    fn test_read_file_format() -> io::Result<()> {
        let data = compress(b"BCF\x02\x01")?;
        let mut reader = Reader::new(&data[..]);

        let (major, minor) = reader.read_file_format()?;

        assert_eq!(major, 2);
        assert_eq!(minor, 1);

        Ok(())
    }

    #[test]
    fn test_read_file_format_with_an_invalid_magic_number() -> io::Result<()> {
        let data = compress(b"BAM\x02\x01")?;
        let mut reader = Reader::new(&data[..]);
        assert!(reader.read_file_format().is_err());
        Ok(())
    }
}
