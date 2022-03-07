use crate::Record;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut Record,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    use crate::reader::record::read_site;

    let l_shared = match reader.read_u32_le().await {
        Ok(n) => usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    let l_indiv = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    buf.resize(l_shared, Default::default());
    reader.read_exact(buf).await?;
    let mut reader = &buf[..];
    let (n_fmt, n_sample) = read_site(&mut reader, record)?;

    let genotypes = record.genotypes_mut().as_mut();
    genotypes.resize(l_indiv, Default::default());
    reader.read_exact(genotypes).await?;
    record.genotypes_mut().set_format_count(n_fmt);
    record.genotypes_mut().set_sample_count(n_sample);

    Ok(l_shared + l_indiv)
}
