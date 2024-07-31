//! BAM record encoder.

mod bin;
mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;

pub(crate) use self::{
    cigar::put_cigar, data::put_data, mapping_quality::put_mapping_quality, name::put_name,
    quality_scores::put_quality_scores, sequence::put_sequence,
};

use std::{error, fmt, io};

use bstr::BStr;
use bytes::BufMut;
use noodles_sam::{
    self as sam,
    alignment::{record_buf::Cigar, Record},
};

use self::{
    bin::put_bin, flags::put_flags, position::put_position,
    reference_sequence_id::put_reference_sequence_id,
};

/// An error when a BAM record fails to encode.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EncodeError {
    /// The reference sequence ID is invalid.
    InvalidReferenceSequenceId(reference_sequence_id::EncodeError),
    /// The alignment start is invalid.
    InvalidAlignmentStart(position::EncodeError),
    /// The mate reference sequence ID is invalid.
    InvalidMateReferenceSequenceId(reference_sequence_id::EncodeError),
    /// The mate alignment start is invalid.
    InvalidMateAlignmentStart(position::EncodeError),
}

impl error::Error for EncodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidReferenceSequenceId(e) => Some(e),
            Self::InvalidAlignmentStart(e) => Some(e),
            Self::InvalidMateReferenceSequenceId(e) => Some(e),
            Self::InvalidMateAlignmentStart(e) => Some(e),
        }
    }
}

impl fmt::Display for EncodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidReferenceSequenceId(_) => write!(f, "invalid reference sequence ID"),
            Self::InvalidAlignmentStart(_) => write!(f, "invalid alignment start"),
            Self::InvalidMateReferenceSequenceId(_) => {
                write!(f, "invalid mate reference sequence ID")
            }
            Self::InvalidMateAlignmentStart(_) => write!(f, "invalid mate alignment start"),
        }
    }
}

pub(crate) fn encode<B, R>(dst: &mut B, header: &sam::Header, record: &R) -> io::Result<()>
where
    B: BufMut,
    R: Record + ?Sized,
{
    // ref_id
    let reference_sequence_id = record.reference_sequence_id(header).transpose()?;
    put_reference_sequence_id(dst, header, reference_sequence_id)
        .map_err(EncodeError::InvalidReferenceSequenceId)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // pos
    let alignment_start = record.alignment_start().transpose()?;
    put_position(dst, alignment_start)
        .map_err(EncodeError::InvalidAlignmentStart)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    put_l_read_name(dst, record.name())?;

    // mapq
    let mapping_quality = record.mapping_quality().transpose()?;
    put_mapping_quality(dst, mapping_quality);

    // bin
    let alignment_end = record.alignment_end().transpose()?;
    put_bin(dst, alignment_start, alignment_end);

    // n_cigar_op
    let cigar = overflowing_put_cigar_op_count(dst, record)?;

    // flag
    let flags = record.flags()?;
    put_flags(dst, flags);

    let l_seq = u32::try_from(record.sequence().len())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    dst.put_u32_le(l_seq);

    // next_ref_id
    let mate_reference_sequence_id = record.mate_reference_sequence_id(header).transpose()?;
    put_reference_sequence_id(dst, header, mate_reference_sequence_id)
        .map_err(EncodeError::InvalidMateReferenceSequenceId)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // next_pos
    let mate_alignment_start = record.mate_alignment_start().transpose()?;
    put_position(dst, mate_alignment_start)
        .map_err(EncodeError::InvalidMateAlignmentStart)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    // tlen
    let template_length = record.template_length()?;
    put_template_length(dst, template_length);

    // read_name
    put_name(dst, record.name())?;

    if let Some(cigar) = &cigar {
        put_cigar(dst, cigar)?;
    } else {
        put_cigar(dst, &record.cigar())?;
    }

    let sequence = record.sequence();
    let base_count = sequence.len();

    // seq
    let read_length = record.cigar().read_length()?;
    put_sequence(dst, read_length, sequence)?;

    // qual
    put_quality_scores(dst, base_count, record.quality_scores())?;

    put_data(dst, record.data())?;

    if cigar.is_some() {
        data::field::put_cigar(dst, &record.cigar())?;
    }

    Ok(())
}

fn put_l_read_name<B>(dst: &mut B, name: Option<&BStr>) -> io::Result<()>
where
    B: BufMut,
{
    use std::mem;

    use self::name::MISSING;

    let mut name_len = name.map(|name| name.len()).unwrap_or(MISSING.len());

    // + NUL terminator
    name_len += mem::size_of::<u8>();

    let l_read_name =
        u8::try_from(name_len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    dst.put_u8(l_read_name);

    Ok(())
}

fn overflowing_put_cigar_op_count<B, R>(dst: &mut B, record: &R) -> io::Result<Option<Cigar>>
where
    B: BufMut,
    R: Record + ?Sized,
{
    use sam::alignment::record::cigar::{op::Kind, Op};

    let cigar = record.cigar();

    if let Ok(op_count) = u16::try_from(cigar.len()) {
        dst.put_u16_le(op_count);
        Ok(None)
    } else {
        dst.put_u16_le(2);

        let k = record.sequence().len();
        let m = cigar.alignment_span()?;

        Ok(Some(
            [Op::new(Kind::SoftClip, k), Op::new(Kind::Skip, m)]
                .into_iter()
                .collect(),
        ))
    }
}

fn put_template_length<B>(dst: &mut B, template_length: i32)
where
    B: BufMut,
{
    dst.put_i32_le(template_length);
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use bstr::BString;
    use noodles_core::Position;
    use noodles_sam::{
        alignment::RecordBuf,
        header::record::value::{map::ReferenceSequence, Map},
    };

    use super::*;

    #[test]
    fn test_encode_with_default_fields() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();
        let header = sam::Header::default();
        let record = RecordBuf::default();
        encode(&mut buf, &header, &record)?;

        let expected = [
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

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_encode_with_all_fields() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::{
            record::{
                cigar::{op::Kind, Op},
                data::field::Tag,
                Flags, MappingQuality,
            },
            record_buf::{data::field::Value, QualityScores, Sequence},
        };

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        const SQ1_LN: NonZeroUsize = match NonZeroUsize::new(13) {
            Some(length) => length,
            None => unreachable!(),
        };

        let mut buf = Vec::new();

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(SQ1_LN))
            .build();

        let record = RecordBuf::builder()
            .set_name(BString::from(b"r0"))
            .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
            .set_reference_sequence_id(1)
            .set_alignment_start(Position::try_from(9)?)
            .set_mapping_quality(MappingQuality::try_from(13)?)
            .set_cigar(
                [Op::new(Kind::Match, 3), Op::new(Kind::SoftClip, 1)]
                    .into_iter()
                    .collect(),
            )
            .set_mate_reference_sequence_id(1)
            .set_mate_alignment_start(Position::try_from(22)?)
            .set_template_length(144)
            .set_sequence(Sequence::from(b"ACGT"))
            .set_quality_scores(QualityScores::from(vec![45, 35, 43, 50]))
            .set_data(
                [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
                    .into_iter()
                    .collect(),
            )
            .build();

        encode(&mut buf, &header, &record)?;

        let expected = [
            0x01, 0x00, 0x00, 0x00, // ref_id = 1
            0x08, 0x00, 0x00, 0x00, // pos = 8
            0x03, // l_read_name = 3
            0x0d, // mapq = 13
            0x49, 0x12, // bin = 4681
            0x02, 0x00, // n_cigar_op = 2
            0x41, 0x00, // flag = 65
            0x04, 0x00, 0x00, 0x00, // l_seq = 4
            0x01, 0x00, 0x00, 0x00, // next_ref_id = 1
            0x15, 0x00, 0x00, 0x00, // next_pos = 21
            0x90, 0x00, 0x00, 0x00, // tlen = 144
            b'r', b'0', 0x00, // read_name = "r0\x00"
            0x30, 0x00, 0x00, 0x00, // cigar[0] = 3M
            0x14, 0x00, 0x00, 0x00, // cigar[1] = 1S
            0x12, 0x48, // seq = ACGT
            0x2d, 0x23, 0x2b, 0x32, // qual = NDLS
            b'N', b'H', b'C', 0x01, // data[0] = NH:i:1
        ];

        assert_eq!(buf, expected);

        Ok(())
    }

    #[test]
    fn test_encode_with_oversized_cigar() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::{
            record::{
                cigar::{op::Kind, Op},
                data::field::Tag,
                Flags,
            },
            record_buf::{data::field::Value, Sequence},
        };

        const BASE_COUNT: usize = 65536;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(131072) {
            Some(length) => length,
            None => unreachable!(),
        };

        let mut buf = Vec::new();

        let header = sam::Header::builder()
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .build();

        let cigar = Cigar::from(vec![Op::new(Kind::Match, 1); BASE_COUNT]);
        let sequence = Sequence::from(vec![b'A'; BASE_COUNT]);

        let record = RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::MIN)
            .set_cigar(cigar)
            .set_sequence(sequence)
            .set_data(
                [(Tag::ALIGNMENT_HIT_COUNT, Value::from(1))]
                    .into_iter()
                    .collect(),
            )
            .build();

        encode(&mut buf, &header, &record)?;

        let mut expected = vec![
            0x00, 0x00, 0x00, 0x00, // ref_id = 0
            0x00, 0x00, 0x00, 0x00, // pos = 1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x49, 0x02, // bin = 585
            0x02, 0x00, // n_cigar_op = 2
            0x00, 0x00, // flag = <empty>
            0x00, 0x00, 0x01, 0x00, // l_seq = 65536
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
            0x04, 0x00, 0x10, 0x00, // cigar[0] = 65536S
            0x03, 0x00, 0x10, 0x00, // cigar[1] = 65536N
        ];

        expected.resize(expected.len() + (BASE_COUNT + 1) / 2, 0x11); // seq = [A, ...]
        expected.resize(expected.len() + BASE_COUNT, 0xff); // qual = [0xff, ...]
        expected.extend([b'N', b'H', b'C', 0x01]); // data[0] = NH:i:1

        // data[1] = CG:B,I:...
        expected.extend([b'C', b'G', b'B', b'I', 0x00, 0x00, 0x01, 0x00]);
        expected.extend((0..BASE_COUNT).flat_map(|_| {
            [
                0x10, 0x00, 0x00, 0x00, // 1M
            ]
        }));

        assert_eq!(buf, expected);

        Ok(())
    }
}
