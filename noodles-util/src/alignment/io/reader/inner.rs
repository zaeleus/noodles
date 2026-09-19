use std::io::{self, BufReader, Read, Seek};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Region;
use noodles_cram as cram;
use noodles_sam as sam;

use crate::alignment::{Index, Record};

pub(super) enum Inner<R> {
    Sam(sam::io::Reader<BufReader<R>>),
    SamGz(sam::io::Reader<bgzf::io::Reader<BufReader<R>>>),
    Bam(bam::io::Reader<bgzf::io::Reader<BufReader<R>>>),
    BamRaw(bam::io::Reader<BufReader<R>>),
    Cram(cram::io::BufReader<BufReader<R>>),
}

impl<R> Inner<R>
where
    R: Read,
{
    pub(super) fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Inner::Sam(reader) => reader.read_header(),
            Inner::SamGz(reader) => reader.read_header(),
            Inner::Bam(reader) => reader.read_header(),
            Inner::BamRaw(reader) => reader.read_header(),
            Inner::Cram(reader) => reader.get_mut().read_header(),
        }
    }

    pub(super) fn read_record(
        &mut self,
        header: &sam::Header,
        record: &mut Record,
    ) -> io::Result<usize> {
        match self {
            Inner::Sam(reader) => {
                if !matches!(record, Record::Sam(_)) {
                    *record = Record::Sam(sam::Record::default());
                }

                if let Record::Sam(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::SamGz(reader) => {
                if !matches!(record, Record::Sam(_)) {
                    *record = Record::Sam(sam::Record::default());
                }

                if let Record::Sam(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::Bam(reader) => {
                if !matches!(record, Record::Bam(_)) {
                    *record = Record::Bam(bam::Record::default());
                }

                if let Record::Bam(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::BamRaw(reader) => {
                if !matches!(record, Record::Bam(_)) {
                    *record = Record::Bam(bam::Record::default());
                }

                if let Record::Bam(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::Cram(reader) => {
                if !matches!(record, Record::Cram(_)) {
                    *record = Record::Cram(sam::alignment::RecordBuf::default());
                }

                if let Record::Cram(r) = record {
                    reader.read_record_buf(header, r)
                } else {
                    unreachable!();
                }
            }
        }
    }

    pub(super) fn records<'r, 'h: 'r>(
        &'r mut self,
        header: &'h sam::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r {
        let records: Box<dyn Iterator<Item = io::Result<_>>> = match self {
            Inner::Sam(reader) => Box::new(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
            Inner::SamGz(reader) => Box::new(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
            Inner::Bam(reader) => Box::new(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
            Inner::BamRaw(reader) => Box::new(reader.records().map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
            Inner::Cram(reader) => Box::new(reader.get_mut().records(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
        };

        records
    }
}

impl<R> Inner<R>
where
    R: Read + Seek,
{
    pub(super) fn query<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i Index,
        region: &Region,
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r> {
        let records: Box<dyn Iterator<Item = io::Result<_>>> = match (self, index) {
            (Inner::SamGz(reader), Index::Sam(idx)) => {
                let query = reader.query(header, idx, region)?;

                Box::new(query.records().map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            (Inner::Bam(reader), Index::Bam(idx)) => {
                let query = reader.query(header, idx, region)?;

                Box::new(query.records().map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            (Inner::Cram(reader), Index::Cram(idx)) => {
                let query = reader.get_mut().query(header, idx, region)?;

                Box::new(query.records().map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "format-index mismatch",
                ));
            }
        };

        Ok(records)
    }

    pub(super) fn query_unmapped<'r, 'h: 'r, 'i: 'r>(
        &'r mut self,
        header: &'h sam::Header,
        index: &'i Index,
    ) -> io::Result<impl Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'r> {
        let records: Box<dyn Iterator<Item = io::Result<_>>> = match (self, index) {
            (Inner::SamGz(reader), Index::Sam(idx)) => {
                let query = reader.query_unmapped(idx)?;

                Box::new(query.map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            (Inner::Bam(reader), Index::Bam(idx)) => {
                let query = reader.query_unmapped(idx)?;

                Box::new(query.map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            (Inner::Cram(reader), Index::Cram(idx)) => {
                let query = reader.get_mut().query_unmapped(header, idx)?;

                Box::new(query.map(|result| {
                    result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
                }))
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "format-index mismatch",
                ));
            }
        };

        Ok(records)
    }
}
