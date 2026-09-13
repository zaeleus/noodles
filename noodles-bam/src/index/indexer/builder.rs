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
        let depth = if let Some(n) = self.max_position_hint {
            fit_depth(n).expect("unsupported max position")
        } else {
            BAI_DEPTH
        };

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

fn fit_depth(max_position: Position) -> Option<u8> {
    let mut depth = BAI_DEPTH;

    loop {
        let max_index_position = calculate_max_position(BAI_MIN_SHIFT, depth)?;

        if max_position <= max_index_position {
            break;
        }

        depth += 1;
    }

    Some(depth)
}

fn calculate_max_position(min_shift: u8, depth: u8) -> Option<Position> {
    const MAX_DEPTH: u8 = 10;

    if depth > MAX_DEPTH {
        None
    } else {
        1u64.checked_shl(u32::from(min_shift) + 3 * u32::from(depth))
            .and_then(|n| usize::try_from(n - 1).ok())
            .and_then(Position::new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fit_depth() {
        #[cfg(not(target_pointer_width = "16"))]
        {
            assert_eq!(fit_depth(Position::MIN), Some(5));

            assert_eq!(
                fit_depth(const { Position::new((1 << 29) - 1).unwrap() }),
                Some(5)
            );
            assert_eq!(
                fit_depth(const { Position::new(1 << 29).unwrap() }),
                Some(6)
            );
            assert_eq!(
                fit_depth(const { Position::new(u32::MAX as usize).unwrap() }),
                Some(6)
            );
        }

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        {
            assert_eq!(
                fit_depth(const { Position::new(1 << 32).unwrap() }),
                Some(7)
            );
            assert_eq!(
                fit_depth(const { Position::new((1 << 44) - 1).unwrap() }),
                Some(10)
            );
            assert!(fit_depth(const { Position::new(1 << 44).unwrap() }).is_none());
            assert!(fit_depth(Position::MAX).is_none());
        }
    }

    #[test]
    fn test_calculate_max_position() {
        #[cfg(not(target_pointer_width = "16"))]
        assert_eq!(
            calculate_max_position(14, 5),
            Some(const { Position::new((1 << 29) - 1).unwrap() })
        );

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        assert_eq!(
            calculate_max_position(14, 10),
            Some(const { Position::new((1 << 44) - 1).unwrap() })
        );

        assert!(calculate_max_position(14, 11).is_none());
    }
}
