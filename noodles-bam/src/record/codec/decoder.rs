//! BAM record decoder.

mod bin;
pub(crate) mod cigar;
pub mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
pub(crate) mod sequence;
mod template_length;

pub(crate) use self::{
    cigar::get_cigar, data::get_data, quality_scores::get_quality_scores,
    reference_sequence_id::get_reference_sequence_id, sequence::get_sequence,
};

use std::{error, fmt};

use bytes::Buf;
use noodles_sam::alignment::RecordBuf;

use self::{
    bin::discard_bin, flags::get_flags, mapping_quality::get_mapping_quality, name::get_name,
    position::get_position, template_length::get_template_length,
};

/// An error when a raw BAM record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DecodeError {
    /// The reference sequence ID is invalid.
    InvalidReferenceSequenceId(reference_sequence_id::DecodeError),
    /// The alignment start is invalid.
    InvalidAlignmentStart(position::DecodeError),
    /// The mapping quality is invalid.
    InvalidMappingQuality(mapping_quality::DecodeError),
    /// The bin is invalid.
    InvalidBin(bin::DecodeError),
    /// The flags are invalid.
    InvalidFlags(flags::DecodeError),
    /// The mate reference sequence ID is invalid.
    InvalidMateReferenceSequenceId(reference_sequence_id::DecodeError),
    /// The mate alignment start is invalid.
    InvalidMateAlignmentStart(position::DecodeError),
    /// The template length is invalid.
    InvalidTemplateLength(template_length::DecodeError),
    /// The name is invalid.
    InvalidName(name::DecodeError),
    /// The CIGAR is invalid.
    InvalidCigar(cigar::DecodeError),
    /// The sequence is invalid.
    InvalidSequence(sequence::DecodeError),
    /// The quality scores are invalid.
    InvalidQualityScores(quality_scores::DecodeError),
    /// The data is invalid.
    InvalidData(data::DecodeError),
}

impl error::Error for DecodeError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidReferenceSequenceId(e) => Some(e),
            Self::InvalidAlignmentStart(e) => Some(e),
            Self::InvalidMappingQuality(e) => Some(e),
            Self::InvalidBin(e) => Some(e),
            Self::InvalidFlags(e) => Some(e),
            Self::InvalidMateReferenceSequenceId(e) => Some(e),
            Self::InvalidMateAlignmentStart(e) => Some(e),
            Self::InvalidTemplateLength(e) => Some(e),
            Self::InvalidName(e) => Some(e),
            Self::InvalidCigar(e) => Some(e),
            Self::InvalidSequence(e) => Some(e),
            Self::InvalidQualityScores(e) => Some(e),
            Self::InvalidData(e) => Some(e),
        }
    }
}

impl fmt::Display for DecodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidReferenceSequenceId(_) => write!(f, "invalid reference sequence ID"),
            Self::InvalidAlignmentStart(_) => write!(f, "invalid alignment start"),
            Self::InvalidMappingQuality(_) => write!(f, "invalid mapping quality"),
            Self::InvalidBin(_) => write!(f, "invalid bin"),
            Self::InvalidFlags(_) => write!(f, "invalid flags"),
            Self::InvalidMateReferenceSequenceId(_) => {
                write!(f, "invalid mate reference sequence ID")
            }
            Self::InvalidMateAlignmentStart(_) => write!(f, "invalid mate alignment start"),
            Self::InvalidTemplateLength(_) => write!(f, "invalid template length"),
            Self::InvalidName(_) => write!(f, "invalid read name"),
            Self::InvalidCigar(_) => write!(f, "invalid CIGAR"),
            Self::InvalidSequence(_) => write!(f, "invalid sequence"),
            Self::InvalidQualityScores(_) => write!(f, "invalid quality scores"),
            Self::InvalidData(_) => write!(f, "invalid data"),
        }
    }
}

pub(crate) fn decode<B>(src: &mut B, record: &mut RecordBuf) -> Result<(), DecodeError>
where
    B: Buf,
{
    *record.reference_sequence_id_mut() =
        get_reference_sequence_id(src).map_err(DecodeError::InvalidReferenceSequenceId)?;

    *record.alignment_start_mut() =
        get_position(src).map_err(DecodeError::InvalidAlignmentStart)?;

    let l_read_name = name::get_length(src).map_err(DecodeError::InvalidName)?;

    *record.mapping_quality_mut() =
        get_mapping_quality(src).map_err(DecodeError::InvalidMappingQuality)?;

    discard_bin(src).map_err(DecodeError::InvalidBin)?;

    let n_cigar_op = cigar::get_op_count(src).map_err(DecodeError::InvalidCigar)?;

    *record.flags_mut() = get_flags(src).map_err(DecodeError::InvalidFlags)?;

    let l_seq = sequence::get_length(src).map_err(DecodeError::InvalidSequence)?;

    *record.mate_reference_sequence_id_mut() =
        get_reference_sequence_id(src).map_err(DecodeError::InvalidMateReferenceSequenceId)?;

    *record.mate_alignment_start_mut() =
        get_position(src).map_err(DecodeError::InvalidMateAlignmentStart)?;

    *record.template_length_mut() =
        get_template_length(src).map_err(DecodeError::InvalidTemplateLength)?;

    get_name(src, record.name_mut(), l_read_name).map_err(DecodeError::InvalidName)?;
    get_cigar(src, record.cigar_mut(), n_cigar_op).map_err(DecodeError::InvalidCigar)?;
    get_sequence(src, record.sequence_mut(), l_seq).map_err(DecodeError::InvalidSequence)?;
    get_quality_scores(src, record.quality_scores_mut(), l_seq)
        .map_err(DecodeError::InvalidQualityScores)?;
    get_data(src, record.data_mut()).map_err(DecodeError::InvalidData)?;

    cigar::resolve(record).map_err(DecodeError::InvalidCigar)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_with_invalid_l_read_name() {
        let data = [
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x00, // l_read_name = 0
        ];
        let mut src = &data[..];

        let mut record = RecordBuf::default();

        assert!(matches!(
            decode(&mut src, &mut record),
            Err(DecodeError::InvalidName(_))
        ));
    }
}
