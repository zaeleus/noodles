use std::io::{self, BufWriter, Write};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Write as _};

pub(super) enum Inner<W>
where
    W: Write,
{
    Bcf(bcf::io::Writer<bgzf::io::Writer<W>>),
    BcfRaw(bcf::io::Writer<BufWriter<W>>),
    Vcf(vcf::io::Writer<BufWriter<W>>),
    VcfGz(vcf::io::Writer<bgzf::io::Writer<W>>),
}

impl<W> Inner<W>
where
    W: Write,
{
    pub(super) fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            Self::Bcf(writer) => writer.write_header(header),
            Self::BcfRaw(writer) => writer.write_header(header),
            Self::Vcf(writer) => writer.write_header(header),
            Self::VcfGz(writer) => writer.write_header(header),
        }
    }

    pub(super) fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        match self {
            Self::Bcf(writer) => writer.write_variant_record(header, record),
            Self::BcfRaw(writer) => writer.write_variant_record(header, record),
            Self::Vcf(writer) => writer.write_variant_record(header, record),
            Self::VcfGz(writer) => writer.write_variant_record(header, record),
        }
    }
}
