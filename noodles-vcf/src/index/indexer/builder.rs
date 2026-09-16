//! VCF indexer builder.

use std::{error, fmt};

use noodles_core::Position;
use noodles_csi::{self as csi, binning_index::index::header::ReferenceSequenceNames};

use super::{Indexer, Inner};
use crate::index::Format;

const TABIX_MIN_SHIFT: u8 = 14;
const TABIX_DEPTH: u8 = 5;

/// An error returned when a VCF indexer fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The estimated max position is unsupported.
    ///
    /// No format is able to support the given estimated max position.
    UnsupportedMaxPosition,
    /// The format is invalid.
    ///
    /// The estimated max position cannot fit in the given format.
    InvalidFormat,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedMaxPosition => write!(f, "unsupported max position"),
            Self::InvalidFormat => write!(f, "invalid format"),
        }
    }
}

/// A VCF indexer builder.
#[derive(Default)]
pub struct Builder {
    format: Option<Format>,
    max_position_hint: Option<Position>,
}

impl Builder {
    /// Sets the VCF index format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::{Format, Indexer};
    /// let builder = Indexer::builder().set_format(Format::Tabix);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets a max position hint.
    ///
    /// This sets the maximum expected position that is used for indexing. When <= 2<sup>29</sup>,
    /// this will select tabix (TBI); otherwise, a coordinate-sorted index (CSI) is used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_core::Position;
    /// use noodles_vcf::index::Indexer;
    ///
    /// let builder = Indexer::builder().set_max_position_hint(Position::MIN);
    /// ```
    pub fn set_max_position_hint(mut self, max_position_hint: Position) -> Self {
        self.max_position_hint = Some(max_position_hint);
        self
    }

    /// Builds a VCF indexer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::Indexer;
    /// let indexer = Indexer::builder().build()?;
    /// # Ok::<_, noodles_vcf::index::indexer::builder::BuildError>(())
    /// ```
    pub fn build(self) -> Result<Indexer, BuildError> {
        let depth = if let Some(n) = self.max_position_hint {
            fit_depth(n).ok_or(BuildError::UnsupportedMaxPosition)?
        } else {
            TABIX_DEPTH
        };

        let format = if depth == TABIX_DEPTH {
            self.format.unwrap_or_default()
        } else if let Some(Format::Tabix) = self.format {
            return Err(BuildError::InvalidFormat);
        } else {
            Format::Csi
        };

        let inner = match format {
            Format::Csi => Inner::Csi(csi::binning_index::Indexer::new(TABIX_MIN_SHIFT, depth)),
            Format::Tabix => Inner::Tabix(csi::binning_index::Indexer::default()),
        };

        Ok(Indexer {
            header: csi::binning_index::index::header::Builder::vcf().build(),
            reference_sequence_names: ReferenceSequenceNames::default(),
            inner,
        })
    }
}

fn fit_depth(max_position: Position) -> Option<u8> {
    let mut depth = TABIX_DEPTH;

    loop {
        let max_index_position = calculate_max_position(TABIX_MIN_SHIFT, depth)?;

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
            .and_then(|n| usize::try_from(n).ok())
            .and_then(Position::new)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build() {
        let builder = Builder::default();
        assert!(builder.build().is_ok());

        #[cfg(not(target_pointer_width = "16"))]
        {
            let builder = Builder::default()
                .set_max_position_hint(const { Position::new((1 << 29) + 1).unwrap() });
            assert!(builder.build().is_ok());

            let builder = Builder::default()
                .set_format(Format::Tabix)
                .set_max_position_hint(const { Position::new((1 << 29) + 1).unwrap() });
            assert_eq!(builder.build().err(), Some(BuildError::InvalidFormat));
        }

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        {
            let builder = Builder::default().set_max_position_hint(Position::MAX);

            assert_eq!(
                builder.build().err(),
                Some(BuildError::UnsupportedMaxPosition)
            );
        }
    }

    #[test]
    fn test_fit_depth() {
        #[cfg(not(target_pointer_width = "16"))]
        {
            assert_eq!(fit_depth(Position::MIN), Some(5));

            assert_eq!(
                fit_depth(const { Position::new(1 << 29).unwrap() }),
                Some(5)
            );
            assert_eq!(
                fit_depth(const { Position::new((1 << 29) + 1).unwrap() }),
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
                Some(6)
            );
            assert_eq!(
                fit_depth(const { Position::new((1 << 32) + 1).unwrap() }),
                Some(7)
            );
            assert_eq!(
                fit_depth(const { Position::new(1 << 44).unwrap() }),
                Some(10)
            );
            assert!(fit_depth(const { Position::new((1 << 44) + 1).unwrap() }).is_none());
            assert!(fit_depth(Position::MAX).is_none());
        }
    }

    #[test]
    fn test_calculate_max_position() {
        #[cfg(not(target_pointer_width = "16"))]
        assert_eq!(
            calculate_max_position(14, 5),
            Some(const { Position::new(1 << 29).unwrap() })
        );

        #[cfg(not(any(target_pointer_width = "16", target_pointer_width = "32")))]
        assert_eq!(
            calculate_max_position(14, 10),
            Some(const { Position::new(1 << 44).unwrap() })
        );

        assert!(calculate_max_position(14, 11).is_none());
    }
}
