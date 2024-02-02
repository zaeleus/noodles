mod bounds;

use std::io;

use self::bounds::Bounds;
use super::Genotypes;

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Fields {
    site_buf: Vec<u8>,
    samples_buf: Vec<u8>,
    bounds: Bounds,
}

impl Fields {
    pub(crate) fn site_buf_mut(&mut self) -> &mut Vec<u8> {
        &mut self.site_buf
    }

    pub(crate) fn samples_buf_mut(&mut self) -> &mut Vec<u8> {
        &mut self.samples_buf
    }

    pub(super) fn reference_sequence_id(&self) -> i32 {
        let src = &self.site_buf[bounds::REFERENCE_SEQUENCE_ID_RANGE];
        // SAFETY: `src` is 4 bytes.
        i32::from_le_bytes(src.try_into().unwrap())
    }

    // N.B. this is 0-based.
    pub(super) fn position(&self) -> i32 {
        let src = &self.site_buf[bounds::POSITION_RANGE];
        // SAFETY: `src` is 4 bytes.
        i32::from_le_bytes(src.try_into().unwrap())
    }

    pub(super) fn span(&self) -> i32 {
        let src = &self.site_buf[bounds::SPAN_RANGE];
        // SAFETY: `src` is 4 bytes.
        i32::from_le_bytes(src.try_into().unwrap())
    }

    pub(super) fn quality_score(&self) -> io::Result<Option<f32>> {
        use crate::record::codec::value::Float;

        let src = &self.site_buf[bounds::QUALITY_SCORE_RANGE];
        // SAFETY: `src` is 4 bytes.
        let n = f32::from_le_bytes(src.try_into().unwrap());

        match Float::from(n) {
            Float::Value(n) => Ok(Some(n)),
            Float::Missing => Ok(None),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid quality score",
            )),
        }
    }

    fn sample_count(&self) -> io::Result<usize> {
        let src = &self.site_buf[bounds::SAMPLE_COUNT_RANGE];
        let n = u32::from_le_bytes([src[0], src[1], src[2], 0x00]);
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    fn format_key_count(&self) -> usize {
        let n = self.site_buf[bounds::FORMAT_KEY_COUNT_INDEX];
        usize::from(n)
    }

    pub(super) fn genotypes(&self) -> io::Result<Genotypes<'_>> {
        self.sample_count().map(|sample_count| {
            Genotypes::new(&self.samples_buf, sample_count, self.format_key_count())
        })
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            site_buf: vec![
                0x00, 0x00, 0x00, 0x00, // chrom = 0
                0x00, 0x00, 0x00, 0x00, // pos = 0 (0-based)
                0x01, 0x00, 0x00, 0x00, // rlen = 1
                0x01, 0x00, 0x80, 0x7f, // qual = None
                0x00, 0x00, // n_info = 0
                0x01, 0x00, // n_allele = 1
                0x00, 0x00, 0x00, // n_sample = 0
                0x00, // n_fmt = 0
                0x07, // ids = []
                0x17, 0x4e, // ref = N
                0x00, // filters = []
            ],
            samples_buf: Vec::new(),
            bounds: Bounds,
        }
    }
}
