use tokio::io::{self, AsyncWrite, AsyncWriteExt};

use super::write_newline;
use crate::record::Sequence;

pub(super) async fn write_sequence<W>(
    writer: &mut W,
    sequence: &Sequence,
    line_base_count: usize,
) -> io::Result<()>
where
    W: AsyncWrite + Unpin,
{
    for bases in sequence.as_ref().chunks(line_base_count) {
        writer.write_all(bases).await?;
        write_newline(writer).await?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_write_sequence() -> io::Result<()> {
        let mut writer = Vec::new();
        let sequence = Sequence::from(b"AC".to_vec());
        write_sequence(&mut writer, &sequence, 4).await?;
        assert_eq!(writer, b"AC\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGT".to_vec());
        write_sequence(&mut writer, &sequence, 4).await?;
        assert_eq!(writer, b"ACGT\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGTACGT".to_vec());
        write_sequence(&mut writer, &sequence, 4).await?;
        assert_eq!(writer, b"ACGT\nACGT\n");

        writer.clear();
        let sequence = Sequence::from(b"ACGTACGTAC".to_vec());
        write_sequence(&mut writer, &sequence, 4).await?;
        assert_eq!(writer, b"ACGT\nACGT\nAC\n");

        Ok(())
    }
}
