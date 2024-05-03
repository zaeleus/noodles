use std::num::NonZeroUsize;

use bstr::BString;
use noodles_sam::{
    self as sam,
    header::{
        record::value::{map::ReferenceSequence, Map},
        ReferenceSequences,
    },
};
use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncRead, AsyncReadExt, BufReader};

use crate::io::reader::{bytes_with_nul_to_bstring, header::reference_sequences_eq};

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let mut header = read_header_inner(reader).await?;
    let reference_sequences = read_reference_sequences(reader).await?;

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

async fn read_header_inner<R>(reader: &mut R) -> io::Result<sam::Header>
where
    R: AsyncRead + Unpin,
{
    let l_text = reader.read_u32_le().await.map(u64::from)?;

    let mut parser = sam::header::Parser::default();

    let mut header_reader = BufReader::new(reader.take(l_text));
    let mut buf = Vec::new();

    while read_line(&mut header_reader, &mut buf).await? != 0 {
        parser
            .parse_partial(&buf)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    discard_padding(&mut header_reader).await?;

    Ok(parser.finish())
}

async fn read_line<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<usize>
where
    R: AsyncBufRead + Unpin,
{
    const NUL: u8 = 0x00;
    const LINE_FEED: u8 = b'\n';
    const CARRIAGE_RETURN: u8 = b'\r';

    dst.clear();

    let src = reader.fill_buf().await?;

    if src.is_empty() || src[0] == NUL {
        return Ok(0);
    }

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

async fn discard_padding<R>(reader: &mut R) -> io::Result<()>
where
    R: AsyncBufRead + Unpin,
{
    loop {
        let src = reader.fill_buf().await?;

        if src.is_empty() {
            return Ok(());
        }

        let len = src.len();
        reader.consume(len);
    }
}

async fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
where
    R: AsyncRead + Unpin,
{
    let n_ref = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref);

    for _ in 0..n_ref {
        let (name, reference_sequence) = read_reference_sequence(reader).await?;
        reference_sequences.insert(name, reference_sequence);
    }

    Ok(reference_sequences)
}

async fn read_reference_sequence<R>(reader: &mut R) -> io::Result<(BString, Map<ReferenceSequence>)>
where
    R: AsyncRead + Unpin,
{
    let l_name = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut c_name = vec![0; l_name];
    reader.read_exact(&mut c_name).await?;

    let name = bytes_with_nul_to_bstring(&c_name)?;

    let l_ref = reader.read_u32_le().await.and_then(|len| {
        usize::try_from(len)
            .and_then(NonZeroUsize::try_from)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let reference_sequence = Map::<ReferenceSequence>::new(l_ref);

    Ok((name, reference_sequence))
}

#[cfg(test)]
mod tests {
    use bytes::BufMut;

    use noodles_sam::header::record::value::{
        map::{self, header::Version},
        Map,
    };

    use super::*;

    const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
        Some(length) => length,
        None => unreachable!(),
    };

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        let mut data = Vec::new();
        data.put_u32_le(27); // l_text
        data.put_slice(b"@HD\tVN:1.6\n@SQ\tSN:sq0\tLN:8\n"); // text
        data.put_u32_le(1); // n_ref
        data.put_u32_le(4); // ref[0].l_name
        data.put_slice(b"sq0\x00"); // ref[0].name
        data.put_u32_le(8); // ref[0].l_ref

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let expected = sam::Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<map::ReferenceSequence>::new(SQ0_LN))
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_reference_sequences() -> io::Result<()> {
        let data = [
            0x01, 0x00, 0x00, 0x00, // n_ref = 1
            0x04, 0x00, 0x00, 0x00, // ref[0].l_name = 4
            0x73, 0x71, 0x30, 0x00, // ref[0].name = "sq0\x00"
            0x08, 0x00, 0x00, 0x00, // ref[0].l_ref = 8
        ];

        let mut reader = &data[..];
        let actual = read_reference_sequences(&mut reader).await?;

        let expected: ReferenceSequences =
            [(BString::from("sq0"), Map::<ReferenceSequence>::new(SQ0_LN))]
                .into_iter()
                .collect();

        assert_eq!(actual, expected);

        Ok(())
    }
}
