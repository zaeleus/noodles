# Changelog

## Unreleased

### Changed

  * csi/async/reader: Disallow duplicate bin IDs.

### Fixed

  * csi/writer: Fix column indices being off-by-one when writing a tabix
    header ([#215]).

[#215]: https://github.com/zaeleus/noodles/issues/215

## 0.26.0 - 2023-10-26

### Changed

  * csi/index/reference_sequence: Change bins to an ordered map ([#213]).

    Bins now maintain their insertion order. While this does not directly
    affect reading and in-memory usage, it does make serialization
    deterministic.

[#213]: https://github.com/zaeleus/noodles/issues/213

## 0.25.1 - 2023-10-19

### Fixed

  * csi/index/indexer: Fix final reference sequence count when building index.

## 0.25.0 - 2023-10-12

### Changed

  * csi: Sync dependencies.

## 0.24.0 - 2023-08-31

### Changed

  * csi: Sync dependencies.

## 0.23.0 - 2023-08-17

### Changed

  * csi: Sync dependencies.

## 0.22.0 - 2023-07-06

### Changed

  * csi: Update to indexmap 2.0.0.

## 0.21.0 - 2023-06-29

### Added

  * csi/index/reference_sequence/bin: Make `Bin::max_id` and `Bin::metadata_id`
    `const`.

### Changed

  * csi/reader: Disallow duplicate bin IDs.

## 0.20.0 - 2023-06-15

### Changed

  * csi: Sync dependencies.

## 0.19.0 - 2023-06-01

### Fixed

  * csi/index/reference_sequence: Use the available linear index for the start
    position of the first record in the last linear bin ([#172]).

[#172]: https://github.com/zaeleus/noodles/issues/172

## 0.18.0 - 2023-05-18

### Added

  * csi/io: Add a filtered indexed records iterator (`FilterByRegion`).

    This filters indexed records that intersect a given region.

  * csi/io: Add an indexed records iterator (`IndexedRecords`).

    This parses lines from a reader as an indexed record.

  * csi/io: Add `IndexedRecord` trait to represent the components used to
    index a record.

    I.e., a reference sequence name, start position, and end position.

  * csi/io: Add an indexed reader (`IndexedReader`).

  * csi/io/query: Add `Query::indexed_records` to create an
    iterator of indexed records.

### Changed

  * csi/index/header: Change column indices to be 0-based.

## 0.17.0 - 2023-05-04

### Added

  * csi/io: Add `Query` reader for reading the uncompressed data between all
    the given chunks.

## 0.16.0 - 2023-04-27

### Changed

  * csi: Sync dependencies.

## 0.15.0 - 2023-04-06

### Added

  * csi/async: Add a convenience write function to write an index to a file
    (`csi::r#async::write`).

  * csi/index: Add an indexer (`csi::index::Indexer`) ([#157]).

  * csi/index: Add `Header`.

    This is the same structure as a tabix header. It was moved from
    `noodles_tabix::index::Header`.

  * csi/index/reference_sequence: Add a linear index
    (`ReferenceSequence::linear_index`).

    This is optional and may be preferred over bin linear offsets.

  * csi/index/reference_sequence/builder: Build linear index.

[#157]: https://github.com/zaeleus/noodles/issues/157

### Changed

  * csi/index: Replace `aux` with an optional tabix header
    (`Option<csi::index::Header>`).

    It isn't clear what `aux` is supposed to be used for, so noodles-csi
    assumes if it's set, it's a tabix header.

  * csi/index: Optimize chunks from `Index::query`.

  * csi/index/reference_sequence: Increase the visibility of `Builder`.

  * csi/index/reference_sequence: Change bins to a `HashMap<usize, Bin>`.

    `Bin` no longer holds its bin ID.

  * csi/index/reference_sequence/bin: Increase the visibility of `Builder`.

### Removed

  * csi: Remove `BinningIndex`.

    `csi::Index` is now the only concrete implementation of a binning index.

## 0.14.0 - 2023-03-03

### Changed

  * csi: Sync dependencies.

## 0.13.0 - 2023-02-03

### Changed

  * csi: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * csi: Remove `BinningIndexReferenceSequence`.

    This was deprecated in noodles-csi 0.4.0. Use
    `noodles_csi::binning_index::ReferenceSequenceExt` instead.

  * csi/index: Remove `Index::unmapped_read_count`.

    This was deprecated in noodles-csi 0.2.0. Use
    `Index::unplaced_unmapped_record_count` instead.

  * csi/index/builder: Remove `Index::set_n_no_coor`.

    This was deprecated in noodles-csi 0.2.0. Use
    `Builder::set_unplaced_unmapped_record_count` instead.

## 0.12.0 - 2022-11-18

### Changed

  * csi: Sync dependencies.

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
