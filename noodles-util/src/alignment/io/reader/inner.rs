use std::io::{self, BufReader, Read};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_sam as sam;

use crate::alignment::Record;

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
