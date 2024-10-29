use noodles_core::Position;

use super::num::write_u16_le;

// § 4.2.1 "BIN field calculation" (2021-06-03): "Note unmapped reads with `POS` 0 (which
// becomes -1 in BAM) therefore use `reg2bin(-1, 0)` which is computed as 4680."
const UNMAPPED_BIN: u16 = 4680;

pub(super) fn write_bin(
    dst: &mut Vec<u8>,
    alignment_start: Option<Position>,
    alignment_end: Option<Position>,
) {
    let n = match (alignment_start, alignment_end) {
        (Some(start), Some(end)) => region_to_bin(start, end),
        _ => UNMAPPED_BIN,
    };

    write_u16_le(dst, n);
}

// § 5.3 "C source code for computing bin number and overlapping bins" (2021-06-03)
#[allow(clippy::eq_op)]
fn region_to_bin(alignment_start: Position, alignment_end: Position) -> u16 {
    let start = usize::from(alignment_start) - 1;
    let end = usize::from(alignment_end) - 1;

    let bin = if start >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (start >> 14)
    } else if start >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (start >> 17)
    } else if start >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (start >> 20)
    } else if start >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (start >> 23)
    } else if start >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (start >> 26)
    } else {
        0
    };

    // SAFETY: Truncate overflowing bin IDs.
    bin as u16
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_bin() {
        let mut buf = Vec::new();

        buf.clear();
        write_bin(&mut buf, None, None);
        assert_eq!(buf, [0x48, 0x12]);

        buf.clear();
        write_bin(&mut buf, Position::new(8), Position::new(13));
        assert_eq!(buf, [0x49, 0x12]);
    }

    #[test]
    fn test_region_to_bin() -> Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;
        assert_eq!(region_to_bin(start, end), 4681);

        let start = Position::try_from(63245986)?;
        let end = Position::try_from(63245986)?;
        assert_eq!(region_to_bin(start, end), 8541);

        // https://github.com/zaeleus/noodles/issues/229
        const DEPTH: usize = 5;
        const MIN_SHIFT: usize = 14;
        const MAX_POSITION: usize = 1 << (MIN_SHIFT + 3 * DEPTH);
        const U16_SIZE: usize = 1 << 16;
        const MAX_BIN_ID: usize = (1 << ((DEPTH + 1) * 3)) / 7;
        const WINDOW_SIZE: usize = 1 << MIN_SHIFT;
        let start =
            Position::try_from((MAX_POSITION + 1) + (WINDOW_SIZE * (U16_SIZE - MAX_BIN_ID)))?;
        let end = start;
        assert_eq!(region_to_bin(start, end), 0);

        Ok(())
    }
}
