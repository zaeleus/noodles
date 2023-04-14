//! BCF reader and iterators.

mod header;
pub(crate) mod query;
pub(crate) mod record;
mod records;
pub(crate) mod string_map;
pub(crate) mod value;

pub use self::{query::Query, records::Records};

use std::io::{self, BufRead, Read, Seek};

use byteorder::ReadBytesExt;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_csi as csi;
use noodles_vcf as vcf;

use self::header::read_header;
use super::Record;
use crate::header::string_maps::{ContigStringMap, StringMaps};

/// A BCF reader.
///
/// The BCF format is comprised of two parts: 1) a VCF header and 2) a list of records.
pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Returns a reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::Reader::from(&data[..]);
    /// assert!(reader.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    /// Returns a mutable reference to the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let mut reader = bcf::Reader::from(&data[..]);
    /// assert!(reader.get_mut().is_empty());
    /// ```
    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    /// Returns the underlying reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::Reader::from(&data[..]);
    /// assert!(reader.into_inner().is_empty());
    /// ```
    pub fn into_inner(self) -> R {
        self.inner
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
        read_magic(&mut self.inner)?;
        read_format_version(&mut self.inner)
    }

    /// Reads the raw VCF header.
    ///
    /// The position of the stream is expected to be directly after the file format.
    ///
    /// This returns the raw VCF header as a [`String`]. It can subsequently be parsed as a
    /// [`noodles_vcf::Header`].
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
        read_header(&mut self.inner)
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
    /// let mut record = bcf::Record::default();
    /// reader.read_record(&mut record)?;
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        record::read_record(&mut self.inner, &mut self.buf, record)
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
    ///     println!("{:?}", record);
    /// }
    /// # Ok::<(), io::Error>(())
    /// ```
    pub fn records(&mut self) -> Records<'_, R> {
        Records::new(self)
    }
}

impl<R> Reader<bgzf::Reader<R>>
where
    R: Read,
{
    /// Creates a BCF reader.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf as bcf;
    /// let data = [];
    /// let reader = bcf::Reader::new(&data[..]);
    /// ```
    pub fn new(reader: R) -> Self {
        Self::from(bgzf::Reader::new(reader))
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

impl<R> Reader<bgzf::Reader<R>>
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
    /// use noodles_bcf::{self as bcf, header::StringMaps};
    /// use noodles_core::Region;
    /// use noodles_csi as csi;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::Reader::new)?;
    /// reader.read_file_format()?;
    ///
    /// let string_maps: StringMaps = reader.read_header()?.parse()?;
    ///
    /// let index = csi::read("sample.bcf.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(string_maps.contigs(), &index, &region)?;
    ///
    /// for result in query {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn query(
        &mut self,
        contig_string_map: &ContigStringMap,
        index: &csi::Index,
        region: &Region,
    ) -> io::Result<Query<'_, R>> {
        let reference_sequence_id = resolve_region(contig_string_map, region)?;
        let chunks = index.query(reference_sequence_id, region.interval())?;

        Ok(Query::new(
            self,
            chunks,
            reference_sequence_id,
            region.interval(),
        ))
    }
}

impl<R> From<R> for Reader<R> {
    fn from(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
    }
}

impl<R> vcf::VariantReader<R> for Reader<R>
where
    R: BufRead,
{
    fn read_variant_header(&mut self) -> io::Result<vcf::Header> {
        read_magic(&mut self.inner)?;
        read_format_version(&mut self.inner)?;

        read_header(&mut self.inner).and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    fn variant_records<'a>(
        &'a mut self,
        header: &'a vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<vcf::Record>> + 'a> {
        let string_maps = StringMaps::from(header);

        Box::new(self.records().map(move |result| {
            result.and_then(|record| record.try_into_vcf_record(header, &string_maps))
        }))
    }
}

fn read_magic<R>(reader: &mut R) -> io::Result<()>
where
    R: Read,
{
    use crate::MAGIC_NUMBER;

    let mut buf = [0; 3];
    reader.read_exact(&mut buf)?;

    if buf == MAGIC_NUMBER {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BCF header",
        ))
    }
}

fn read_format_version<R>(reader: &mut R) -> io::Result<(u8, u8)>
where
    R: Read,
{
    let major_version = reader.read_u8()?;
    let minor_version = reader.read_u8()?;

    Ok((major_version, minor_version))
}

pub(crate) fn resolve_region(
    contig_string_map: &ContigStringMap,
    region: &Region,
) -> io::Result<usize> {
    contig_string_map
        .get_index_of(region.name())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("region does not exist in contigs: {region:?}"),
            )
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_magic() {
        let data = b"BCF";
        let mut reader = &data[..];
        assert!(read_magic(&mut reader).is_ok());

        let data = [];
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        let data = b"BAM";
        let mut reader = &data[..];
        assert!(matches!(
            read_magic(&mut reader),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_read_format_version() -> io::Result<()> {
        let data = [0x02, 0x01];
        let mut reader = &data[..];
        assert_eq!(read_format_version(&mut reader)?, (2, 1));
        Ok(())
    }
}
