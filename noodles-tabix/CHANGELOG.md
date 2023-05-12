# Changelog

## Unreleased

### Changed

  * tabix/async/reader: Ensure header column indices are > 0.

  * tabix/reader: Ensure header column indices are > 0.

## 0.20.0 - 2023-05-04

### Changed

  * tabix: Sync dependencies.

## 0.19.0 - 2023-04-27

### Changed

  * tabix: Sync dependencies.

## 0.18.0 - 2023-04-06

### Changed

  * tabix: Replace `Index` with `noodles_csi::Index`.

### Removed

  * tabix/index: Remove `Header`.

    Use `noodles_csi::index::Header` instead.

## 0.17.0 - 2023-03-03

### Changed

  * tabix: Sync dependencies.

## 0.16.0 - 2023-02-03

### Added

  * tabix/index/header/format: Implement `std::error::Error::source` for
    `TryFromIntError`.

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * tabix/index: Remove `ReferenceSequenceNames`.

    This was deprecated in noodles-tabix 0.11.0. Use
    `noodles_tabix::index::header::ReferenceSequenceNames` instead.

  * tabix/index: Remove `Index::reference_sequence_names`.

    This was deprecated in noodles-tabix 0.11.0. Use
    `Header::reference_sequence_names` instead.

  * tabix/index: Remove `Index::unmapped_read_count`.

    This was deprecated in noodles-tabix 0.2.0. Use
    `Index::unplaced_unmapped_record_count` instead.

  * tabix/index/builder: Remove `Builder::set_unmapped_read_count`.

    This was deprecated in noodles-tabix 0.3.0. Use
    `Builder::set_unplaced_unmapped_record_count` instead.

## 0.15.0 - 2022-11-18

### Changed

  * tabix: Sync dependencies.

## 0.14.0 - 2022-10-28

### Changed

  * tabix: Sync dependencies.

## 0.13.0 - 2022-10-20

### Changed

  * tabix: Sync dependencies.

## 0.12.1 - 2022-09-29

### Fixed

  * tabix: Sync dependencies.

## 0.12.0 - 2022-08-16

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.11.0 - 2022-07-05

### Changed

  * tabix/index/reference_sequence/bin: Change bin ID to a `usize`.

  * tabix/index: Move reference sequence names to `Header`.

### Deprecated

  * tabix/index: Deprecate `ReferenceSequenceNames`.

    Use `header::ReferenceSequenceNames` instead.

  * tabix/index: Deprecate `Index::reference_sequence_names`.

    Use `Header::reference_sequence_names` instead.

## 0.10.0 - 2022-06-09

### Changed

  * tabix/reader: Fail if reference sequence names buffer has trailing data.

### Fixed

  * tabix: Sync dependencies.

## 0.9.1 - 2022-06-08

### Fixed

  * tabix: Sync dependencies.

## 0.9.0 - 2022-03-29

### Changed

  * tabix/index/indexer: Change start and end to `Position`.

  * tabix/index/reference_sequence: `ReferenceSequence::query` returns an
    `io::Error` instead of `QueryError`.

### Fixed

  * tabix/index/reference_sequence: Ensure the start position is not out of
    range for a query (`2^29 - 1`).

## 0.8.1 - 2022-03-02

### Fixed

  * tabix: Sync dependencies.

## 0.8.0 - 2022-02-17

### Added

  * tabix: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.7.3 - 2022-01-27

### Fixed

  * tabix: Sync dependencies.

## 0.7.2 - 2021-12-02

### Fixed

  * tabix: Require tokio's `fs` feature as a dependency ([#62]).

[#62]: https://github.com/zaeleus/noodles/issues/62

## 0.7.1 - 2021-11-18

### Fixed

  * tabix: Sync dependencies.

## 0.7.0 - 2021-11-11

### Changed

  * tabix: Update to Rust 2021.

## 0.6.1 - 2021-09-19

### Fixed

  * tabix/async: Fix writer not finalizing.

## 0.6.0 - 2021-08-19

### Changed

  * tabix: Update to tokio 1.10.0.

### Fixed

  * tabix: Define features to enable for Docs.rs.

## 0.5.0 - 2021-08-11

### Added

  * tabix: Add convenience write function to write an index to a file:
    `tabix::write`.

  * tabix/async: Add async reader (`tabix::AsyncReader`).

  * tabix/async: Add async writer (`tabix::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.4.0 - 2021-08-04

### Changed

  * tabix/reader: Disallow duplicate reference sequence names.

## 0.3.0 - 2021-07-30

### Deprecated

  * tabix/index/builder: Deprecate `set_unmapped_read_count`.

    Use `set_unplaced_unmapped_record_count` instead.

### Fixed

  * tabix/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * tabix/index: Implemented `BinningIndex` for `Index`.

  * tabix/index: Added `query` method to find chunks that intersect the given
    region.

  * tabix/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

### Changed

  * tabix/index: Reference sequence names are stored as an
    `index::ReferenceSequenceNames` (`IndexSet<String>`).

### Deprecated

  * tabix/index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

### Fixed

  * tabix: Fixed documentation link in package manifest ([#31]).

  * tabix/reader: Avoid casts that may truncate.

    Fields that convert from `i32` to other integer types now check whether
    they are in range.

  * tabix/writer: Avoid casts that may truncate.

    Fields that convert to `i32` from other integer types now check whether
    they are in range.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * tabix/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * tabix: Initial release.
