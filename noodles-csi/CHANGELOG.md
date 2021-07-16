# Changelog

## Unreleased

### Added

  * Add convenience function to write an entire index to a file: `csi::write`.

  * binning_index: Added chunk merging functions for chunk list reduction
    (`noodles_csi::binning_index::{merge_chunks, optimize_chunks}`).

    Chunks are merged when they overlap and can be filtered by a minimum
    offset.

  * binning_index: Added `BinningIndex` and `BinningIndexReferenceSequence`
    traits to define shared behavior among binning index formats.

## 0.1.0 - 2021-07-14

  * Initial release.
