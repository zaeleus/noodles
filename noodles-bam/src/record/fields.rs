//! BAM record fields.

mod bounds;

use std::{io, mem};

use bstr::{BStr, ByteSlice};

use self::bounds::Bounds;
use super::{Cigar, Data, QualityScores, Sequence};

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds,
}

impl Fields {
    pub(super) fn reference_sequence_id(&self) -> Option<i32> {
        let src = &self.buf[bounds::REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
    }

    // N.B. this is 0-based.
    pub(super) fn alignment_start(&self) -> Option<i32> {
        let src = &self.buf[bounds::ALIGNMENT_START_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_position(src.try_into().unwrap())
    }

    pub(super) fn mapping_quality(&self) -> Option<u8> {
        const MISSING: u8 = 255;

        match self.buf[bounds::MAPPING_QUALITY_INDEX] {
            MISSING => None,
            n => Some(n),
        }
    }

    pub(super) fn flags(&self) -> u16 {
        let src = &self.buf[bounds::FLAGS_RANGE];
        // SAFETY: `src` is 2 bytes.
        u16::from_le_bytes(src.try_into().unwrap())
    }

    pub(super) fn mate_reference_sequence_id(&self) -> Option<i32> {
        let src = &self.buf[bounds::MATE_REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        get_reference_sequence_id(src.try_into().unwrap())
    }

    pub(super) fn mate_alignment_start(&self) -> Option<i32> {
        let src = &self.buf[bounds::MATE_ALIGNMENT_START_RANGE];
        get_position(src.try_into().unwrap())
    }

    pub(super) fn template_length(&self) -> i32 {
        let src = &self.buf[bounds::TEMPLATE_LENGTH_RANGE];
        // SAFETY: `src` is 4 bytes.
        i32::from_le_bytes(src.try_into().unwrap())
    }

    pub(super) fn name(&self) -> Option<&BStr> {
        const NUL: u8 = 0x00;
        const MISSING: &[u8] = &[b'*', NUL];

        match &self.buf[self.bounds.name_range()] {
            MISSING => None,
            buf => Some(buf.strip_suffix(&[NUL]).unwrap_or(buf).as_bstr()),
        }
    }

    pub(super) fn cigar(&self) -> Cigar<'_> {
        use super::data::get_raw_cigar;

        const SKIP: u8 = 3;
        const SOFT_CLIP: u8 = 4;

        fn decode_op(buf: &[u8]) -> (u8, usize) {
            // SAFETY: `buf` is 4 bytes.
            let n = u32::from_le_bytes(buf.try_into().unwrap());
            ((n & 0x0f) as u8, usize::try_from(n >> 4).unwrap())
        }

        let src = &self.buf[self.bounds.cigar_range()];

        if src.len() == 2 * mem::size_of::<u32>() {
            let k = self.sequence().len();

            // SAFETY: `src` is 8 bytes.
            let op_1 = decode_op(&src[0..4]);
            let op_2 = decode_op(&src[4..8]);

            if op_1 == (SOFT_CLIP, k) && matches!(op_2, (SKIP, _)) {
                let mut data_src = &self.buf[self.bounds.data_range()];

                if let Ok(Some(buf)) = get_raw_cigar(&mut data_src) {
                    return Cigar::new(buf);
                }
            }
        }

        Cigar::new(src)
    }

    pub(super) fn sequence(&self) -> Sequence<'_> {
        let src = &self.buf[self.bounds.sequence_range()];

        let buf = &self.buf[bounds::READ_LENGTH_RANGE];
        // SAFETY: `buf.len() == 4`.
        let n = u32::from_le_bytes(buf.try_into().unwrap());
        let base_count = usize::try_from(n).unwrap();

        Sequence::new(src, base_count)
    }

    pub(super) fn quality_scores(&self) -> QualityScores<'_> {
        const MISSING: u8 = 0xff;

        let buf = &self.buf[self.bounds.quality_scores_range()];

        // § 4.2.3 "SEQ and QUAL encoding" (2024-11-06): "When base quality are omitted but the
        // sequence is not, `qual` is filled with `0xFF` bytes (to length `l_seq`)."
        let src = if buf.iter().all(|&b| b == MISSING) {
            &[]
        } else {
            buf
        };

        QualityScores::new(src)
    }

    pub(super) fn data(&self) -> Data<'_> {
        let src = &self.buf[self.bounds.data_range()];
        Data::new(src)
    }

    pub(crate) fn index(&mut self) -> io::Result<()> {
        index(&self.buf[..], &mut self.bounds)
    }
}

impl Default for Fields {
    fn default() -> Self {
        let buf = vec![
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
            b'*', 0x00, // read_name = "*\x00"
        ];

        let bounds = Bounds {
            name_end: buf.len(),
            cigar_end: buf.len(),
            sequence_end: buf.len(),
            quality_scores_end: buf.len(),
        };

        Self { buf, bounds }
    }
}

impl TryFrom<Vec<u8>> for Fields {
    type Error = io::Error;

    fn try_from(buf: Vec<u8>) -> Result<Self, Self::Error> {
        let mut fields = Self {
            buf,
            bounds: Bounds {
                name_end: 0,
                cigar_end: 0,
                sequence_end: 0,
                quality_scores_end: 0,
            },
        };

        fields.index()?;

        Ok(fields)
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

fn index(buf: &[u8], bounds: &mut Bounds) -> io::Result<()> {
    const MIN_BUF_LENGTH: usize = bounds::TEMPLATE_LENGTH_RANGE.end;

    if buf.len() < MIN_BUF_LENGTH {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let read_name_len = usize::from(buf[bounds::NAME_LENGTH_INDEX]);
    bounds.name_end = bounds::TEMPLATE_LENGTH_RANGE.end + read_name_len;

    let src = &buf[bounds::CIGAR_OP_COUNT_RANGE];
    // SAFETY: `src` is 2 bytes.
    let cigar_op_count = usize::from(u16::from_le_bytes(src.try_into().unwrap()));
    let cigar_len = mem::size_of::<u32>() * cigar_op_count;
    bounds.cigar_end = bounds.name_end + cigar_len;

    let src = &buf[bounds::READ_LENGTH_RANGE];
    // SAFETY: `src` is 4 bytes.
    let base_count = usize::try_from(u32::from_le_bytes(src.try_into().unwrap()))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let sequence_len = base_count.div_ceil(2);
    bounds.sequence_end = bounds.cigar_end + sequence_len;

    bounds.quality_scores_end = bounds.sequence_end + base_count;

    if buf.len() < bounds.quality_scores_end {
        Err(io::Error::from(io::ErrorKind::UnexpectedEof))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static DATA: &[u8] = &[
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

    #[test]
    fn test_name() -> io::Result<()> {
        let fields = Fields::try_from(Vec::from(DATA))?;
        assert!(fields.name().is_none());
        Ok(())
    }

    #[test]
    fn test_name_with_name() -> io::Result<()> {
        let data = vec![
            0xff, 0xff, 0xff, 0xff, // ref_id = -1
            0xff, 0xff, 0xff, 0xff, // pos = -1
            0x02, // l_read_name = 3
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

        let fields = Fields::try_from(data)?;
        assert_eq!(fields.name(), Some(b"r0".as_bstr()));

        Ok(())
    }

    #[test]
    fn test_cigar() -> io::Result<()> {
        let fields = Fields::try_from(Vec::from(DATA))?;
        let cigar = fields.cigar();
        assert_eq!(cigar.as_ref(), &DATA[34..38]);
        Ok(())
    }

    #[test]
    fn test_cigar_with_2_cigar_ops() -> io::Result<()> {
        let data = [
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

        let fields = Fields::try_from(Vec::from(&data))?;
        let cigar = fields.cigar();
        assert_eq!(cigar.as_ref(), &data[34..42]);

        Ok(())
    }

    #[test]
    fn test_cigar_with_overflowing_cigar() -> io::Result<()> {
        let data = [
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

        let fields = Fields::try_from(Vec::from(&data))?;
        let cigar = fields.cigar();
        assert_eq!(cigar.as_ref(), &data[56..]);

        Ok(())
    }

    #[test]
    fn test_index() -> io::Result<()> {
        let mut fields = Fields::default();

        fields.buf.clear();
        fields.buf.extend(DATA);

        fields.index()?;

        assert_eq!(fields.bounds.name_range(), 32..34);
        assert_eq!(fields.bounds.cigar_range(), 34..38);
        assert_eq!(fields.bounds.sequence_range(), 38..40);
        assert_eq!(fields.bounds.quality_scores_range(), 40..44);
        assert_eq!(fields.bounds.data_range(), 44..);

        Ok(())
    }
}
