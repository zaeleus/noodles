//! Async BAM record field readers.

pub mod data;

use std::mem;

use noodles_sam as sam;
use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::{
    record::{Cigar, Data, QualityScores, Sequence},
    Record,
};

pub(super) async fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let block_size = match reader.read_u32_le().await {
        Ok(bs) => usize::try_from(bs).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    record.ref_id = reader.read_i32_le().await?;
    record.pos = reader.read_i32_le().await?;

    let l_read_name = reader.read_u8().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    record.mapq = read_mapping_quality(reader).await?;

    record.bin = reader.read_u16_le().await?;

    let n_cigar_op = reader.read_u16_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    *record.flags_mut() = read_flag(reader).await?;

    let l_seq = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    record.next_ref_id = reader.read_i32_le().await?;
    record.next_pos = reader.read_i32_le().await?;
    record.tlen = reader.read_i32_le().await?;

    read_read_name(reader, &mut record.read_name, l_read_name).await?;
    read_cigar(reader, record.cigar_mut(), n_cigar_op).await?;
    read_seq(reader, record.sequence_mut(), l_seq).await?;
    read_qual(reader, record.quality_scores_mut(), l_seq).await?;
    read_data(
        reader,
        record.data_mut(),
        block_size,
        l_read_name,
        n_cigar_op,
        l_seq,
    )
    .await?;

    Ok(block_size)
}

async fn read_mapping_quality<R>(reader: &mut R) -> io::Result<sam::record::MappingQuality>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_u8()
        .await
        .map(sam::record::MappingQuality::from)
}

async fn read_flag<R>(reader: &mut R) -> io::Result<sam::record::Flags>
where
    R: AsyncRead + Unpin,
{
    reader.read_u16_le().await.map(sam::record::Flags::from)
}

async fn read_read_name<R>(
    reader: &mut R,
    read_name: &mut Vec<u8>,
    l_read_name: usize,
) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    read_name.resize(l_read_name, Default::default());
    reader.read_exact(read_name).await?;
    Ok(())
}

async fn read_cigar<R>(reader: &mut R, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    cigar.resize(n_cigar_op, Default::default());
    cigar.clear();

    for _ in 0..n_cigar_op {
        let op = reader.read_u32_le().await?;
        cigar.push(op);
    }

    Ok(())
}

async fn read_seq<R>(reader: &mut R, seq: &mut Sequence, l_seq: usize) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    seq.set_base_count(l_seq);

    let seq_len = (l_seq + 1) / 2;
    seq.resize(seq_len, Default::default());
    reader.read_exact(seq).await?;

    Ok(())
}

async fn read_qual<R>(reader: &mut R, qual: &mut QualityScores, l_seq: usize) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    qual.resize(l_seq, Default::default());
    reader.read_exact(qual).await?;
    Ok(())
}

async fn read_data<R>(
    reader: &mut R,
    data: &mut Data,
    block_size: usize,
    l_read_name: usize,
    n_cigar_op: usize,
    l_seq: usize,
) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let cigar_len = mem::size_of::<u32>() * n_cigar_op;
    let seq_len = (l_seq + 1) / 2;
    let data_offset = 32 + l_read_name + cigar_len + seq_len + l_seq;
    let data_len = block_size - data_offset;

    data.resize(data_len, Default::default());
    reader.read_exact(data).await?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_record() -> io::Result<()> {
        let data = [
            0x22, 0x00, 0x00, 0x00, // block_size = 34
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x00, 0x00, // n_cigar_op = 0
            0x04, 0x00, // flag = 4
            0x00, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            0x2a, 0x00, // read_name = "*\x00"
        ];

        let mut reader = &data[..];
        let mut record = Record::default();
        let block_size = read_record(&mut reader, &mut record).await?;

        assert_eq!(block_size, 34);
        assert_eq!(record, Record::default());

        Ok(())
    }
}
