use std::io::{self, BufReader, Read};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_sam as sam;

pub(super) enum Inner<R> {
    Sam(sam::io::Reader<BufReader<R>>),
    SamGz(sam::io::Reader<bgzf::io::Reader<BufReader<R>>>),
    Bam(bam::io::Reader<bgzf::io::Reader<BufReader<R>>>),
    BamRaw(bam::io::Reader<BufReader<R>>),
    Cram(cram::io::Reader<BufReader<R>>),
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
            Inner::Cram(reader) => reader.read_header(),
        }
    }

    pub(super) fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn sam::alignment::Record>>> + 'a {
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
            Inner::Cram(reader) => Box::new(reader.records(header).map(|result| {
                result.map(|record| Box::new(record) as Box<dyn sam::alignment::Record>)
            })),
        };

        records
    }
}
