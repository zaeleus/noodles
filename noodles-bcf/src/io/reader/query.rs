mod records;

use std::io;

use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
use noodles_vcf::{self as vcf, variant::Record as _};

use self::records::Records;
use super::Reader;
use crate::Record;

/// A reader over records of a BCF reader that intersects a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h, R> {
    reader: Reader<csi::io::Query<'r, R>>,
    header: &'h vcf::Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<'r, 'h, R> Query<'r, 'h, R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub(super) fn new(
        reader: &'r mut R,
        header: &'h vcf::Header,
        chunks: Vec<Chunk>,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader: Reader::from(csi::io::Query::new(reader, chunks)),
            header,
            reference_sequence_id,
            interval,
        }
    }

    /// Reads a record.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    /// use noodles_csi as csi;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let index = csi::fs::read("sample.bcf.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// let mut record = bcf::Record::default();
    ///
    /// while reader.read_record(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        next_record(
            &mut self.reader,
            record,
            self.header,
            self.reference_sequence_id,
            self.interval,
        )
    }

    /// Returns an iterator over records.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::fs::File;
    /// use noodles_bcf as bcf;
    /// use noodles_csi as csi;
    ///
    /// let mut reader = File::open("sample.bcf").map(bcf::io::Reader::new)?;
    /// let header = reader.read_header()?;
    ///
    /// let index = csi::fs::read("sample.bcf.csi")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn records(self) -> Records<'r, 'h, R> {
        Records::new(self)
    }
}

fn intersects(
    header: &vcf::Header,
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let reference_sequence_name = record.reference_sequence_name(header.string_maps())?;

    let id = header
        .string_maps()
        .contigs()
        .get_index_of(reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "reference sequence name does not exist in contigs: {reference_sequence_name}"
                ),
            )
        })?;

    if reference_sequence_id != id {
        return Ok(false);
    }

    if interval_is_unbounded(region_interval) {
        Ok(true)
    } else {
        let Some(start) = record.variant_start().transpose()? else {
            return Ok(false);
        };

        let end = record.variant_end(header)?;
        let record_interval = Interval::from(start..=end);

        Ok(record_interval.intersects(region_interval))
    }
}

fn interval_is_unbounded(interval: Interval) -> bool {
    interval.start().is_none() && interval.end().is_none()
}

fn next_record<R>(
    reader: &mut Reader<csi::io::Query<'_, R>>,
    record: &mut Record,
    header: &vcf::Header,
    reference_sequence_id: usize,
    interval: Interval,
) -> io::Result<usize>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    loop {
        match reader.read_record(record)? {
            0 => return Ok(0),
            n => {
                if intersects(header, record, reference_sequence_id, interval)? {
                    return Ok(n);
                }
            }
        }
    }
}
