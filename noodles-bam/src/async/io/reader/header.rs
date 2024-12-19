mod magic_number;
mod reference_sequences;
mod sam_header;

use noodles_sam::{self as sam, header::ReferenceSequences};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt};

use self::{magic_number::read_magic_number, reference_sequences::read_reference_sequences};
use crate::{io::reader::header::reference_sequences_eq, MAGIC_NUMBER};

struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: AsyncRead + Unpin,
{
    fn new(inner: R) -> Self {
        Self { inner }
    }

    async fn read_magic_number(&mut self) -> io::Result<[u8; MAGIC_NUMBER.len()]> {
        read_magic_number(&mut self.inner).await
    }

    async fn raw_sam_header_reader(&mut self) -> io::Result<sam_header::Reader<&mut R>> {
        let len = self.inner.read_u32_le().await.map(u64::from)?;
        Ok(sam_header::Reader::new(&mut self.inner, len))
    }

    async fn read_reference_sequences(&mut self) -> io::Result<ReferenceSequences> {
        read_reference_sequences(&mut self.inner).await
    }
}

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let mut header_reader = Reader::new(reader);
    read_header_inner(&mut header_reader).await
}

async fn read_header_inner<R>(reader: &mut Reader<R>) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::header::magic_number;

    reader
        .read_magic_number()
        .await
        .and_then(magic_number::validate)?;

    let mut raw_sam_header_reader = reader.raw_sam_header_reader().await?;
    let mut header = read_sam_header(&mut raw_sam_header_reader).await?;

    let reference_sequences = reader.read_reference_sequences().await?;

    if header.reference_sequences().is_empty() {
        *header.reference_sequences_mut() = reference_sequences;
    } else if !reference_sequences_eq(header.reference_sequences(), &reference_sequences) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "SAM header and binary reference sequence dictionaries mismatch",
        ));
    }

    Ok(header)
}

async fn read_sam_header<R>(reader: &mut sam_header::Reader<R>) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let mut parser = sam::header::Parser::default();

    let mut buf = Vec::new();

    while read_line(reader, &mut buf).await? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    reader.discard_to_end().await?;

    Ok(parser.finish())
}

async fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    match reader.read_until(LINE_FEED, dst).await? {
        0 => Ok(0),
        n => {
            if dst.ends_with(&[LINE_FEED]) {
                dst.pop();

                if dst.ends_with(&[CARRIAGE_RETURN]) {
                    dst.pop();
                }
            }

            Ok(n)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use noodles_sam::header::record::value::{
        map::{self, header::Version},
        Map,
    };

    use super::*;
    use crate::MAGIC_NUMBER;

    const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
        Some(length) => length,
        None => unreachable!(),
    };

    fn put_u32_le(buf: &mut Vec<u8>, n: u32) {
        buf.extend(n.to_le_bytes());
    }

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        let mut src = Vec::new();
        src.extend(MAGIC_NUMBER);
        put_u32_le(&mut src, 27); // l_text
        src.extend(b"@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n"); // text
        put_u32_le(&mut src, 1); // n_ref
        put_u32_le(&mut src, 4); // ref[0].l_name
        src.extend(b"sq0\x00"); // ref[0].name
        put_u32_le(&mut src, 8); // ref[0].l_ref

        let mut reader = &src[..];
        let actual = read_header(&mut reader).await?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
