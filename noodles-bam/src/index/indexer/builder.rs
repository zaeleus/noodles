use noodles_core::Position;
use noodles_csi as csi;

use super::{Indexer, Inner};
use crate::index::Format;

const BAI_MIN_SHIFT: u8 = 14;
const BAI_DEPTH: u8 = 5;

/// A BAM indexer builder.
#[derive(Default)]
pub struct Builder {
    format: Format,
    max_position_hint: Option<Position>,
}

impl Builder {
    /// Sets the BAM index format.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Sets a max position hint.
    ///
    /// This optionally sets the maximum expected position that is used for indexing. When < 2^29,
    /// this will select a BAM index (BAI); otherwise, a coordinate-sorted index (CSI) is used.
    pub fn set_max_position_hint(mut self, max_position_hint: Position) -> Self {
        self.max_position_hint = Some(max_position_hint);
        self
    }

    /// Builds a BAM indexer.
    pub fn build(self) -> Indexer {
        let depth = self.max_position_hint.map(fit_depth).unwrap_or(BAI_DEPTH);

        let format = if depth == BAI_DEPTH {
            self.format
        } else {
            Format::Csi
        };

        let inner = match format {
            Format::Bai => Inner::Bai(csi::binning_index::Indexer::default()),
            Format::Csi => Inner::Csi(csi::binning_index::Indexer::new(BAI_MIN_SHIFT, depth)),
        };

        Indexer { inner }
    }
}

fn fit_depth(max_position: Position) -> u8 {
    const MAX_DEPTH: u8 = 10;

    let mut depth = BAI_DEPTH;

    loop {
        let max_index_position = calculate_max_position(BAI_MIN_SHIFT, depth);

        if max_position <= max_index_position || depth >= MAX_DEPTH {
            break;
        }

        depth += 1;
    }

    depth
}

fn calculate_max_position(min_shift: u8, depth: u8) -> Position {
    let n = (1 << (usize::from(min_shift) + 3 * usize::from(depth))) - 1;
    Position::new(n).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_depth() {
        #[cfg(not(target_pointer_width = "16"))]
        {
            assert_eq!(fit_depth(Position::MIN), 5);

            assert_eq!(
                fit_depth(const { Position::new((1 << 29) - 1).unwrap() }),
                5
            );
            assert_eq!(fit_depth(const { Position::new(1 << 29).unwrap() }), 6);
            assert_eq!(
                fit_depth(const { Position::new(u32::MAX as usize).unwrap() }),
                6
            );
        }

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        {
            assert_eq!(fit_depth(const { Position::new(1 << 32).unwrap() }), 7);
            assert_eq!(fit_depth(Position::MAX), 10);
        }
    }

    #[test]
    fn test_calculate_max_position() {
        assert_eq!(
            calculate_max_position(14, 5),
            const { Position::new(536870911).unwrap() }
        )
    }
}
