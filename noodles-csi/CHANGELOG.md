# Changelog

## Unreleased

### Changed

  * Update to Rust 2021.

### Deprecated

  * csi/binning_index: Rename `csi::BinningIndexReferenceSequence` to
    `csi::binning_index::ReferenceSequenceExt`.

## 0.3.0 - 2021-08-19

### Added

  * async: Add async reader (`csi::AsyncReader`).

  * async: Add async writer (`csi::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.2 - 2021-08-11

### Fixed

  * Sync dependencies.

## 0.2.1 - 2021-07-30

### Fixed

  * reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

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
