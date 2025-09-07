mod reference_sequence;

use noodles_sam::header::ReferenceSequences;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use self::reference_sequence::read_reference_sequence;

pub(super) async fn read_reference_sequences<R>(reader: &mut R) -> io::Result<ReferenceSequences>
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

#[cfg(test)]
mod tests {
    use std::num::NonZero;

    use bstr::BString;
    use noodles_sam::header::record::value::{Map, map::ReferenceSequence};

    use super::*;

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

        let expected: ReferenceSequences = [(
            BString::from("sq0"),
            Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
        )]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        Ok(())
    }
}
