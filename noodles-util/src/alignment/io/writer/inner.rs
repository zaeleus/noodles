use std::io::{self, BufWriter, Write};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_sam::{self as sam, alignment::io::Write as _};

pub(super) enum Inner<W>
where
    W: Write,
{
    Sam(sam::io::Writer<BufWriter<W>>),
    SamGz(sam::io::Writer<bgzf::io::Writer<W>>),
    Bam(bam::io::Writer<bgzf::io::Writer<W>>),
    BamRaw(bam::io::Writer<BufWriter<W>>),
    Cram(cram::io::Writer<W>),
}

impl<W> Inner<W>
where
    W: Write,
{
    pub(super) fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.write_header(header),
            Self::SamGz(writer) => writer.write_header(header),
            Self::Bam(writer) => writer.write_header(header),
            Self::BamRaw(writer) => writer.write_header(header),
            Self::Cram(writer) => writer.write_header(header),
        }
    }

    pub(super) fn write_record<R>(&mut self, header: &sam::Header, record: &R) -> io::Result<()>
    where
        R: sam::alignment::Record,
    {
        match self {
            Self::Sam(writer) => writer.write_alignment_record(header, record),
            Self::SamGz(writer) => writer.write_alignment_record(header, record),
            Self::Bam(writer) => writer.write_alignment_record(header, record),
            Self::BamRaw(writer) => writer.write_alignment_record(header, record),
            Self::Cram(writer) => writer.write_alignment_record(header, record),
        }
    }

    pub(super) fn finish(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.finish(header),
            Self::SamGz(writer) => writer.finish(header),
            Self::Bam(writer) => writer.finish(header),
            Self::BamRaw(writer) => writer.finish(header),
            Self::Cram(writer) => writer.finish(header),
        }
    }
}
