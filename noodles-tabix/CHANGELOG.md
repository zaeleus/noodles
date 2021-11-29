# Changelog

## Unreleased

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
