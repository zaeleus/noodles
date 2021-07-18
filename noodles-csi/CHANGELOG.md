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

  * binning_index: Added `first_record_in_last_linear_bin_start_position`.

    This is the closest position to the unplaced, unmapped records, if any,
    that is available in an index.

  * index: Implemented `BinningIndex` for `Index`.

  * index: Added `query` method to find chunks that intersect the given region.

  * index/reference_sequence: Implemented `BinningIndexReferenceSequence` for
    `ReferenceSequence`.

### Deprecated

  * index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

  * index/builder: Deprecated `Builder::set_n_no_coor`.

    Use `set_unplaced_unmapped_record_count` instead.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
