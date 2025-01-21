use noodles_csi::binning_index::index::Header;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_header<R>(reader: &mut R) -> io::Result<Header>
where
    R: AsyncRead + Unpin,
{
    use std::mem;

    use noodles_csi::io::reader::index::read_header;

    const STATIC_FIELD_COUNT: usize = 6;
    const STATIC_HEADER_SIZE: usize = mem::size_of::<i32>() * STATIC_FIELD_COUNT;

    let mut buf = vec![0; STATIC_HEADER_SIZE];
    reader.read_exact(&mut buf).await?;

    let l_nm = reader.read_i32_le().await?;

    let names_len =
        usize::try_from(l_nm).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut names = vec![0; names_len];
    reader.read_exact(&mut names).await?;

    buf.extend(l_nm.to_le_bytes());
    buf.extend(names);

    let mut buf_reader = &buf[..];
    read_header(&mut buf_reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use bstr::BString;
    use noodles_csi as csi;

    use super::*;

    #[tokio::test]
    async fn test_read_header() -> io::Result<()> {
        let data = [
            0x00, 0x00, 0x00, 0x00, // format = Generic(GFF)
            0x01, 0x00, 0x00, 0x00, // col_seq = 1
            0x04, 0x00, 0x00, 0x00, // col_beg = 4
            0x05, 0x00, 0x00, 0x00, // col_end = 5
            0x23, 0x00, 0x00, 0x00, // meta = '#'
            0x00, 0x00, 0x00, 0x00, // skip = 0
            0x08, 0x00, 0x00, 0x00, // l_nm = 8
            b's', b'q', b'0', 0x00, // names[0] = "sq0"
            b's', b'q', b'1', 0x00, // names[1] = "sq1"
        ];

        let mut reader = &data[..];
        let actual = read_header(&mut reader).await?;

        let names = [BString::from("sq0"), BString::from("sq1")]
            .into_iter()
            .collect();
        let expected = csi::binning_index::index::header::Builder::gff()
            .set_reference_sequence_names(names)
            .build();

        assert_eq!(actual, expected);

        Ok(())
    }
}
