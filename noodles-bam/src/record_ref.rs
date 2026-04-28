use std::{io, mem, ops::Range};

use bstr::{BStr, ByteSlice};
use noodles_core::Position;
use noodles_sam::alignment::record::{Flags, MappingQuality};

use super::record::{
    Cigar, Data, QualityScores, Sequence, try_to_position, try_to_reference_sequence_id,
};

const ALIGNMENT_START_RANGE: Range<usize> = 4..8;
const NAME_LENGTH_INDEX: usize = 8;
const MAPPING_QUALITY_INDEX: usize = 9;
const CIGAR_OP_COUNT_RANGE: Range<usize> = 12..14;
const FLAGS_RANGE: Range<usize> = 14..16;
const READ_LENGTH_RANGE: Range<usize> = 16..20;
const MATE_REFERENCE_SEQUENCE_ID_RANGE: Range<usize> = 20..24;
const MATE_ALIGNMENT_START_RANGE: Range<usize> = 24..28;
const TEMPLATE_LENGTH_RANGE: Range<usize> = 28..32;

const HEAD_SIZE: usize = TEMPLATE_LENGTH_RANGE.end;

/// An immutable view over a BAM record.
pub struct RecordRef<'a> {
    head: &'a [u8; HEAD_SIZE],
    rest: &'a [u8],
}

impl<'a> RecordRef<'a> {
    /// Creates an immutable view over a a BAM record.
    ///
    /// The input must be at minimum 32 bytes; otherwise, `None` is returned. It is the
    /// responsibility of the caller to ensure the input is record-like.
    pub fn new(src: &'a [u8]) -> Option<Self> {
        src.split_first_chunk()
            .map(|(head, rest)| Self { head, rest })
    }

    pub(crate) fn new_unchecked(src: &'a [u8]) -> Self {
        let (head, rest) = src.split_at(HEAD_SIZE);

        Self {
            head: head.try_into().unwrap(),
            rest,
        }
    }

    /// Returns the reference sequence ID.
    pub fn reference_sequence_id(&self) -> Option<io::Result<usize>> {
        // SAFETY: `self.head.len() >= mem::size_of::<i32>()`.
        let src = self.head.first_chunk().unwrap();
        get_reference_sequence_id(*src).map(try_to_reference_sequence_id)
    }

    /// Returns the alignment start.
    pub fn alignment_start(&self) -> Option<io::Result<Position>> {
        let src = &self.head[ALIGNMENT_START_RANGE];
        // SAFETY: `src.len() == mem::size_of::<i32>()`.
        get_position(src.try_into().unwrap()).map(try_to_position)
    }

    fn name_length(&self) -> usize {
        let n = &self.head[NAME_LENGTH_INDEX];
        usize::from(*n)
    }

    /// Returns the mapping quality.
    pub fn mapping_quality(&self) -> Option<MappingQuality> {
        let n = self.head[MAPPING_QUALITY_INDEX];
        MappingQuality::new(n)
    }

    fn cigar_op_count(&self) -> usize {
        let src = &self.head[CIGAR_OP_COUNT_RANGE];
        // SAFETY: `src.len() == mem::size_of::<u16>()`.
        usize::from(u16::from_le_bytes(src.try_into().unwrap()))
    }

    /// Returns the flags.
    pub fn flags(&self) -> Flags {
        let src = &self.head[FLAGS_RANGE];
        // SAFETY: `src.len() == mem::size_of::<u16>()`.
        let n = u16::from_le_bytes(src.try_into().unwrap());
        Flags::from(n)
    }

    pub(crate) fn base_count(&self) -> usize {
        let src = &self.head[READ_LENGTH_RANGE];
        // SAFETY: `src.len() == mem::size_of::<u32>()`.
        let n = u32::from_le_bytes(src.try_into().unwrap());
        usize::try_from(n).unwrap()
    }

    /// Returns the mate reference sequence ID.
    pub fn mate_reference_sequence_id(&self) -> Option<io::Result<usize>> {
        let src = &self.head[MATE_REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src.len() == mem::size_of::<i32>()`.
        get_reference_sequence_id(src.try_into().unwrap()).map(try_to_reference_sequence_id)
    }

    /// Returns the mate alignment start.
    pub fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        let src = &self.head[MATE_ALIGNMENT_START_RANGE];
        // SAFETY: `src.len() == mem::size_of::<i32>()`.
        get_position(src.try_into().unwrap()).map(try_to_position)
    }

    /// Returns the template length.
    pub fn template_length(&self) -> i32 {
        let src = &self.head[TEMPLATE_LENGTH_RANGE];
        // SAFETY: `src.len() == mem::size_of::<i32>()`.
        i32::from_le_bytes(src.try_into().unwrap())
    }

    /// Returns the read name.
    pub fn name(&self) -> Option<&'a BStr> {
        const NUL: u8 = 0x00;
        const MISSING: &[u8] = &[b'*', NUL];

        let end = self.name_length();

        match &self.rest[..end] {
            MISSING => None,
            buf => Some(buf.strip_suffix(&[NUL]).unwrap_or(buf).as_bstr()),
        }
    }

    /// Returns the CIGAR operations.
    pub fn cigar(&self) -> Cigar<'a> {
        use crate::record::data::get_raw_cigar;

        const SKIP: u8 = 3;
        const SOFT_CLIP: u8 = 4;

        fn decode_op(buf: &[u8; 4]) -> (u8, usize) {
            let n = u32::from_le_bytes(*buf);
            ((n & 0x0f) as u8, usize::try_from(n >> 4).unwrap())
        }

        let start = self.name_length();
        let end = start + (self.cigar_op_count() * mem::size_of::<u32>());
        let src = &self.rest[start..end];

        if let ([chunk_0, chunk_1], []) = src.as_chunks() {
            let k = self.base_count();

            let op_1 = decode_op(chunk_0);
            let op_2 = decode_op(chunk_1);

            if op_1 == (SOFT_CLIP, k) && matches!(op_2, (SKIP, _)) {
                let mut data_src = self.raw_data();

                if let Ok(Some(buf)) = get_raw_cigar(&mut data_src) {
                    return Cigar::new(buf);
                }
            }
        }

        Cigar::new(src)
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> Sequence<'a> {
        let (src, base_count) = self.raw_sequence();
        Sequence::new(src, base_count)
    }

    fn raw_sequence(&self) -> (&'a [u8], usize) {
        let start = self.name_length() + (self.cigar_op_count() * mem::size_of::<u32>());

        let base_count = self.base_count();
        let sequence_len = base_count.div_ceil(2);
        let end = start + sequence_len;

        (&self.rest[start..end], base_count)
    }

    /// Returns the quality scores.
    pub fn quality_scores(&self) -> QualityScores<'a> {
        QualityScores::new(self.raw_quality_scores())
    }

    fn raw_quality_scores(&self) -> &'a [u8] {
        const MISSING: u8 = 0xff;

        let base_count = self.base_count();

        let start = self.name_length()
            + (self.cigar_op_count() * mem::size_of::<u32>())
            + base_count.div_ceil(2);

        let end = start + base_count;

        let src = &self.rest[start..end];

        // § 4.2.3 "SEQ and QUAL encoding" (2024-11-06): "When base quality are omitted but the
        // sequence is not, `qual` is filled with `0xFF` bytes (to length `l_seq`)."
        if src.iter().all(|&b| b == MISSING) {
            &[]
        } else {
            src
        }
    }

    /// Returns the data.
    pub fn data(&self) -> Data<'a> {
        Data::new(self.raw_data())
    }

    fn raw_data(&self) -> &'a [u8] {
        let base_count = self.base_count();

        let start = self.name_length()
            + (self.cigar_op_count() * mem::size_of::<u32>())
            + base_count.div_ceil(2)
            + base_count;

        &self.rest[start..]
    }
}

fn get_reference_sequence_id(src: [u8; 4]) -> Option<i32> {
    const UNMAPPED: i32 = -1;

    match i32::from_le_bytes(src) {
        UNMAPPED => None,
        n => Some(n),
    }
}

fn get_position(src: [u8; 4]) -> Option<i32> {
    const MISSING: i32 = -1;

    match i32::from_le_bytes(src) {
        MISSING => None,
        n => Some(n),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fields() -> io::Result<()> {
        const SRC: &[u8; 44] = &[
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x01, 0x00, // n_cigar_op = 1
            0x04, 0x00, // flag = 4
            0x04, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
            0x40, 0x00, 0x00, 0x00, // cigar = 4M
            0x12, 0x48, // sequence = ACGT
            b'N', b'D', b'L', b'S', // quality scores
        ];

        let record = RecordRef::new_unchecked(SRC);

        assert!(record.reference_sequence_id().transpose()?.is_none());
        assert!(record.alignment_start().transpose()?.is_none());
        assert!(record.mapping_quality().is_none());
        assert_eq!(record.flags(), Flags::UNMAPPED);
        assert!(record.mate_reference_sequence_id().transpose()?.is_none());
        assert!(record.mate_alignment_start().transpose()?.is_none());
        assert_eq!(record.template_length(), 0);
        assert!(record.name().is_none());
        assert_eq!(record.cigar().as_ref(), [0x40, 0x00, 0x00, 0x00]);
        assert_eq!(record.sequence().as_ref(), &[0x12, 0x48]);
        assert_eq!(record.quality_scores().as_ref(), b"NDLS");
        assert!(record.data().is_empty());

        Ok(())
    }

    #[test]
    fn test_name() -> io::Result<()> {
        const SRC: &[u8; 45] = &[
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x03, // l_read_name = 3
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x01, 0x00, // n_cigar_op = 1
            0x04, 0x00, // flag = 4
            0x04, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'r', b'0', 0x00, // read_name = "r0\x00"
            0x40, 0x00, 0x00, 0x00, // cigar = 4M
            0x12, 0x48, // sequence = ACGT
            b'N', b'D', b'L', b'S', // quality scores
        ];

        let record = RecordRef::new_unchecked(SRC);
        assert_eq!(record.name(), Some(b"r0".as_bstr()));

        Ok(())
    }

    #[test]
    fn test_cigar_with_2_cigar_ops() -> io::Result<()> {
        const SRC: &[u8; 48] = &[
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x02, 0x00, // n_cigar_op = 2
            0x04, 0x00, // flag = 4
            0x04, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
            0x20, 0x00, 0x00, 0x00, 0x20, 0x00, 0x00, 0x00, // cigar = 2M2M
            0x12, 0x48, // sequence = ACGT
            b'N', b'D', b'L', b'S', // quality scores
        ];

        let record = RecordRef::new_unchecked(SRC);

        assert_eq!(
            record.cigar().as_ref(),
            [0x20, 0x00, 0x00, 0x00, 0x20, 0x00, 0x00, 0x00]
        );

        Ok(())
    }

    #[test]
    fn test_cigar_with_overflowing_cigar() -> io::Result<()> {
        const SRC: &[u8; 60] = &[
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x02, 0x00, // n_cigar_op = 2
            0x04, 0x00, // flag = 4
            0x04, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
            0x44, 0x00, 0x00, 0x00, 0x23, 0x00, 0x00, 0x00, // cigar = 2S2N
            0x12, 0x48, // sequence = ACGT
            b'N', b'D', b'L', b'S', // quality scores
            b'C', b'G', b'B', b'I', 0x01, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00,
            0x00, // data["CG"] = [4M]
        ];

        let record = RecordRef::new_unchecked(SRC);
        assert_eq!(record.cigar().as_ref(), [0x40, 0x00, 0x00, 0x00]);

        Ok(())
    }

    #[test]
    fn test_quality_scores_with_missing_scores() -> io::Result<()> {
        const SRC: &[u8; 44] = &[
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 2
            0xff, // mapq = 255
            0x48, 0x12, // bin = 4680
            0x01, 0x00, // n_cigar_op = 1
            0x04, 0x00, // flag = 4
            0x04, 0x00, 0x00, 0x00, // l_seq = 0
            0xff, 0xff, 0xff, 0xff, // next_ref_id = -1
            0xff, 0xff, 0xff, 0xff, // next_pos = -1
            0x00, 0x00, 0x00, 0x00, // tlen = 0
            b'*', 0x00, // read_name = "*\x00"
            0x40, 0x00, 0x00, 0x00, // cigar = 4M
            0x12, 0x48, // sequence = ACGT
            0xff, 0xff, 0xff, 0xff, // quality scores
        ];

        let record = RecordRef::new_unchecked(SRC);
        assert!(record.quality_scores().is_empty());

        Ok(())
    }
}
