use std::io;

use bytes::BufMut;
use noodles_sam as sam;

use super::record::{
    put_bin, put_cigar, put_data, put_flags, put_l_read_name, put_mapping_quality, put_position,
    put_quality_scores, put_read_name, put_sequence, put_template_length,
};

// § 4.2.3 SEQ and QUAL encoding (2021-06-03)
pub(crate) const NULL_QUALITY_SCORE: u8 = 255;

pub fn encode_alignment_record<B>(
    dst: &mut B,
    header: &sam::Header,
    record: &dyn sam::AlignmentRecord,
) -> io::Result<()>
where
    B: BufMut,
{
    // ref_id
    let reference_sequence_name = record
        .reference_sequence(header.reference_sequences())
        .transpose()?
        .map(|rs| rs.name());

    put_reference_sequence_id(dst, header.reference_sequences(), reference_sequence_name)?;

    // pos
    put_position(dst, record.alignment_start())?;

    put_l_read_name(dst, record.read_name())?;

    // mapq
    put_mapping_quality(dst, record.mapping_quality());

    // bin
    put_bin(dst, record.alignment_start(), record.alignment_end())?;

    let n_cigar_op = u16::try_from(record.cigar().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u16_le(n_cigar_op);

    // flag
    put_flags(dst, record.flags());

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u32_le(l_seq);

    // next_ref_id
    let mate_reference_sequence_name = record
        .mate_reference_sequence(header.reference_sequences())
        .transpose()?
        .map(|rs| rs.name());

    put_reference_sequence_id(
        dst,
        header.reference_sequences(),
        mate_reference_sequence_name,
    )?;

    // next_pos
    put_position(dst, record.mate_alignment_start())?;

    // tlen
    put_template_length(dst, record.template_length());

    put_read_name(dst, record.read_name());

    put_cigar(dst, record.cigar())?;

    // § 4.2.3 SEQ and QUAL encoding (2021-06-03)
    let sequence = record.sequence();
    let quality_scores = record.quality_scores();

    // seq
    put_sequence(dst, sequence);

    if sequence.len() == quality_scores.len() {
        put_quality_scores(dst, quality_scores);
    } else if quality_scores.is_empty() {
        for _ in 0..sequence.len() {
            dst.put_u8(NULL_QUALITY_SCORE);
        }
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "quality scores length mismatch: expected {}, got {}",
                sequence.len(),
                quality_scores.len()
            ),
        ));
    }

    put_data(dst, record.data())?;

    Ok(())
}

fn put_reference_sequence_id<B>(
    dst: &mut B,
    reference_sequences: &sam::header::ReferenceSequences,
    reference_sequence_name: Option<&sam::record::ReferenceSequenceName>,
) -> io::Result<()>
where
    B: BufMut,
{
    use crate::record::reference_sequence_id;

    let id = match reference_sequence_name {
        Some(name) => reference_sequences
            .get_index_of(name.as_str())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence name: {}", name),
                )
            })
            .and_then(|i| {
                i32::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?,
        None => reference_sequence_id::UNMAPPED,
    };

    dst.put_i32_le(id);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        use sam::header::{reference_sequence, ReferenceSequence};

        let reference_sequences = [("sq0".parse()?, 8), ("sq1".parse()?, 13)]
            .into_iter()
            .map(|(name, len): (reference_sequence::Name, i32)| {
                let sn = name.to_string();
                ReferenceSequence::new(name, len).map(|rs| (sn, rs))
            })
            .collect::<Result<_, _>>()?;

        let mut buf = Vec::new();

        buf.clear();
        let reference_sequence_name = "sq0".parse()?;
        put_reference_sequence_id(
            &mut buf,
            &reference_sequences,
            Some(&reference_sequence_name),
        )?;
        assert_eq!(buf, [0x00, 0x00, 0x00, 0x00]);

        buf.clear();
        put_reference_sequence_id(&mut buf, &reference_sequences, None)?;
        assert_eq!(buf, [0xff, 0xff, 0xff, 0xff]);

        buf.clear();
        let reference_sequence_name = "sq2".parse()?;
        assert!(matches!(
            put_reference_sequence_id(&mut buf, &reference_sequences, Some(&reference_sequence_name)),
            Err(ref e) if e.kind() == io::ErrorKind::InvalidInput,
        ));

        Ok(())
    }
}
