use std::io::{self, Read};

use noodles_vcf::{self as vcf, variant::RecordBuf};

use crate::io::reader::num::read_u32_le;

pub(super) fn read_record_buf<R>(
    reader: &mut R,
    header: &vcf::Header,
    buf: &mut Vec<u8>,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: Read,
{
    use crate::record::codec::decoder::{read_samples, read_site};

    let l_shared = match read_u32_le(reader) {
        Ok(n) => usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    let l_indiv = read_u32_le(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    buf.resize(l_shared, 0);
    reader.read_exact(buf)?;
    let mut src = &buf[..];
    let (n_fmt, n_sample) = read_site(&mut src, header, record)?;

    buf.resize(l_indiv, 0);
    reader.read_exact(buf)?;
    let mut src = &buf[..];

    *record.samples_mut() = read_samples(&mut src, header, n_sample, n_fmt)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok(l_shared + l_indiv)
}
