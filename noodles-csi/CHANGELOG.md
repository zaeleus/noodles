# Changelog

## 0.11.0 - 2022-10-28

### Changed

  * csi: Sync dependencies.

## 0.10.0 - 2022-10-20

### Changed

  * csi: Sync dependencies.

## 0.9.1 - 2022-09-29

### Fixed

  * csi: Sync dependencies.

## 0.9.0 - 2022-08-16

### Changed

  * csi: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.8.0 - 2022-07-05

### Added

  * csi/index/reference_sequence: Add finding the minimum start virtual
    position for a bin (`ReferenceSequence::min_offset`).

  * csi/index/reference_sequence/bin/chunk: Implement
    `From<Range<bgzf::VirtualPosition>>` for `Chunk`.

### Fixed

  * csi/index/reference_sequence: Fix finding the start position of the first
    record in the last linear bin.

    The bins are not necessarily sorted.

## 0.7.0 - 2022-06-08

### Changed

  * csi/binning_index: Change generic reference sequence type to an associated
    type.

  * csi/binning_index: Change query interval to `Into<Interval>`.

  * csi/index/reference_sequence/bin: Change bin ID to a `usize`.

## 0.6.0 - 2022-03-29

### Changed

  * csi/binning_index: Change interval bounds to `Position`.

  * csi/index: Change `min_shift` and `depth` to a `u8`.

  * csi/index/reference_sequence: `ReferenceSequence::query` returns an
    `io::Error` instead of `QueryError`.

### Fixed

  * csi/index/reference_sequence: Ensure the start position is not out of range
    for a query.

## 0.5.1 - 2022-03-02

### Fixed

  * csi: Sync dependencies.

## 0.5.0 - 2022-02-17

### Added

  * csi: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.4.3 - 2022-01-27

### Fixed

  * csi: Sync dependencies.

## 0.4.2 - 2021-12-02

### Fixed

  * csi: Require tokio's `fs` feature as a dependency ([#62]).

[#62]: https://github.com/zaeleus/noodles/issues/62

## 0.4.1 - 2021-11-18

### Fixed

  * csi: Sync dependencies.

## 0.4.0 - 2021-11-11

### Changed

  * csi: Update to Rust 2021.

### Deprecated

  * csi/binning_index: Rename `csi::BinningIndexReferenceSequence` to
    `csi::binning_index::ReferenceSequenceExt`.

## 0.3.0 - 2021-08-19

### Added

  * csi/async: Add async reader (`csi::AsyncReader`).

  * csi/async: Add async writer (`csi::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.2 - 2021-08-11

### Fixed

  * csi: Sync dependencies.

## 0.2.1 - 2021-07-30

### Fixed

  * csi/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * csi: Add convenience function to write an entire index to a file:
    `csi::write`.

  * csi/binning_index: Added chunk merging functions for chunk list reduction
    (`noodles_csi::binning_index::{merge_chunks, optimize_chunks}`).

    Chunks are merged when they overlap and can be filtered by a minimum
    offset.

  * csi/binning_index: Added `BinningIndex` and `BinningIndexReferenceSequence`
    traits to define shared behavior among binning index formats.

  * csi/binning_index: Added `first_record_in_last_linear_bin_start_position`.

    This is the closest position to the unplaced, unmapped records, if any,
    that is available in an index.

  * csi/index: Implemented `BinningIndex` for `Index`.

  * csi/index: Added `query` method to find chunks that intersect the given
    region.

  * csi/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

### Deprecated

  * csi/index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

  * csi/index/builder: Deprecated `Builder::set_n_no_coor`.

    Use `set_unplaced_unmapped_record_count` instead.

### Fixed

  * csi: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * csi: Initial release.
