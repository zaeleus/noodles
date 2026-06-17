mod records;

use std::{
    io::{self, Read, Seek, SeekFrom},
    slice, vec,
};

use noodles_core::region::Interval;
use noodles_sam as sam;

use self::records::Records;
use super::{Container, Reader};
use crate::crai;

/// A reader over records that intersect a given region.
///
/// This is created by calling [`Reader::query`].
pub struct Query<'r, 'h: 'r, 'i: 'r, R>
where
    R: Read + Seek,
{
    reader: &'r mut Reader<R>,

    header: &'h sam::Header,

    index: slice::Iter<'i, crai::Record>,

    reference_sequence_id: usize,
    interval: Interval,

    records: vec::IntoIter<sam::alignment::RecordBuf>,
}

impl<'r, 'h: 'r, 'i: 'r, R> Query<'r, 'h, 'i, R>
where
    R: Read + Seek,
{
    pub(super) fn new(
        reader: &'r mut Reader<R>,
        header: &'h sam::Header,
        index: &'i crai::Index,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> Self {
        Self {
            reader,

            header,

            index: index.iter(),

            reference_sequence_id,
            interval,

            records: Vec::new().into_iter(),
        }
    }

    /// Reads a record into an alignment record buffer.
    ///
    /// If successful, 1 is returned. If 0 is returned, the stream reached EOF.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    /// use noodles_fasta as fasta;
    /// use noodles_sam as sam;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let mut query = reader.query(&header, &index, &region)?;
    ///
    /// let mut record = sam::alignment::RecordBuf::default();
    ///
    /// while query.read_record_buf(&mut record)? != 0 {
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn read_record_buf(&mut self, record: &mut sam::alignment::RecordBuf) -> io::Result<usize> {
        loop {
            match self.records.next() {
                Some(r) => {
                    if intersects(&r, self.interval) {
                        *record = r;
                        return Ok(1);
                    }
                }
                None => {
                    if self.read_next_container().transpose()?.is_none() {
                        return Ok(0);
                    }
                }
            }
        }
    }

    /// Returns an iterator over records.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::{fs::File, io};
    /// use noodles_cram::{self as cram, crai};
    /// use noodles_fasta as fasta;
    ///
    /// let mut reader = File::open("sample.cram").map(cram::io::Reader::new)?;
    ///
    /// let header = reader.read_header()?;
    /// let index = crai::fs::read("sample.cram.crai")?;
    /// let region = "sq0:8-13".parse()?;
    /// let query = reader.query(&header, &index, &region)?;
    ///
    /// for result in query.records() {
    ///     let record = result?;
    ///     // ...
    /// }
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn records(self) -> Records<'r, 'h, 'i, R> {
        Records::new(self)
    }

    fn read_next_container(&mut self) -> Option<io::Result<()>> {
        let index_record = self.index.next()?;

        if index_record.reference_sequence_id() != Some(self.reference_sequence_id) {
            return Some(Ok(()));
        }

        if let Err(e) = self.reader.seek(SeekFrom::Start(index_record.offset())) {
            return Some(Err(e));
        }

        let mut container = Container::default();

        match self.reader.read_container(&mut container) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(e)),
        };

        let compression_header = match container.compression_header() {
            Ok(compression_header) => compression_header,
            Err(e) => return Some(Err(e)),
        };

        let records = container
            .slices()
            .map(|result| {
                let slice = result?;

                let (core_data_src, external_data_srcs) = slice.decode_blocks()?;

                slice
                    .records(
                        self.reader.reference_sequence_repository.clone(),
                        self.header,
                        &compression_header,
                        &core_data_src,
                        &external_data_srcs,
                    )
                    .and_then(|records| {
                        records
                            .into_iter()
                            .map(|record| {
                                sam::alignment::RecordBuf::try_from_alignment_record(
                                    self.header,
                                    &record,
                                )
                            })
                            .collect::<io::Result<Vec<_>>>()
                    })
            })
            .collect::<Result<Vec<_>, _>>();

        let records = match records {
            Ok(records) => records,
            Err(e) => return Some(Err(e)),
        };

        self.records = records
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .into_iter();

        Some(Ok(()))
    }
}

fn intersects(record: &sam::alignment::RecordBuf, region_interval: Interval) -> bool {
    match (record.alignment_start(), record.alignment_end()) {
        (Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}
