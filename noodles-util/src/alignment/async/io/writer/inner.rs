use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_sam as sam;
use tokio::io::{self, AsyncWrite, AsyncWriteExt, BufWriter};

pub(super) enum Inner<W>
where
    W: AsyncWrite,
{
    Sam(sam::r#async::io::Writer<BufWriter<W>>),
    SamGz(sam::r#async::io::Writer<bgzf::r#async::io::Writer<W>>),
    Bam(bam::r#async::io::Writer<bgzf::r#async::io::Writer<W>>),
    BamRaw(bam::r#async::io::Writer<BufWriter<W>>),
    Cram(cram::r#async::io::Writer<W>),
}

impl<W> Inner<W>
where
    W: AsyncWrite + Unpin,
{
    pub(super) async fn write_header(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.write_header(header).await,
            Self::SamGz(writer) => writer.write_header(header).await,
            Self::Bam(writer) => writer.write_header(header).await,
            Self::BamRaw(writer) => writer.write_header(header).await,
            Self::Cram(writer) => writer.write_header(header).await,
        }
    }

    pub(super) async fn write_record(
        &mut self,
        header: &sam::Header,
        record: &dyn sam::alignment::Record,
    ) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.write_alignment_record(header, record).await,
            Self::SamGz(writer) => writer.write_alignment_record(header, record).await,
            Self::Bam(writer) => writer.write_alignment_record(header, record).await,
            Self::BamRaw(writer) => writer.write_alignment_record(header, record).await,
            Self::Cram(writer) => writer.write_alignment_record(header, record).await,
        }
    }

    pub(super) async fn shutdown(&mut self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Sam(writer) => writer.get_mut().shutdown().await,
            Self::SamGz(writer) => writer.get_mut().shutdown().await,
            Self::Bam(writer) => writer.shutdown().await,
            Self::BamRaw(writer) => writer.shutdown().await,
            Self::Cram(writer) => writer.shutdown(header).await,
        }
    }
}
