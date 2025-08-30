use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;
use tokio::io::{self, AsyncWrite, BufWriter};

pub(super) enum Inner<W>
where
    W: AsyncWrite,
{
    Bcf(bcf::r#async::io::Writer<bgzf::r#async::io::Writer<W>>),
    BcfRaw(bcf::r#async::io::Writer<BufWriter<W>>),
    Vcf(vcf::r#async::io::Writer<BufWriter<W>>),
    VcfGz(vcf::r#async::io::Writer<bgzf::r#async::io::Writer<W>>),
}

impl<W> Inner<W>
where
    W: AsyncWrite + Unpin,
{
    pub(super) async fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        match self {
            Self::Bcf(writer) => writer.write_header(header).await,
            Self::BcfRaw(writer) => writer.write_header(header).await,
            Self::Vcf(writer) => writer.write_header(header).await,
            Self::VcfGz(writer) => writer.write_header(header).await,
        }
    }

    pub async fn write_record(
        &mut self,
        header: &vcf::Header,
        record: &dyn vcf::variant::Record,
    ) -> io::Result<()> {
        match self {
            Self::Bcf(writer) => writer.write_variant_record(header, record).await,
            Self::BcfRaw(writer) => writer.write_variant_record(header, record).await,
            Self::Vcf(writer) => writer.write_variant_record(header, record).await,
            Self::VcfGz(writer) => writer.write_variant_record(header, record).await,
        }
    }
}
