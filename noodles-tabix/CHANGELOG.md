# Changelog

## Unreleased

### Added

  * index: Implemented `BinningIndex` for `Index`.

  * index: Added `query` method to find chunks that intersect the given region.

  * index/reference_sequence: Implemented `BinningIndexReferenceSequence` for
    `ReferenceSequence`.

### Changed

  * index: Reference sequence names are stored as an
    `index::ReferenceSequenceNames` (`IndexSet<String>`).

### Deprecated

  * index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

  * reader: Avoid casts that may truncate.

    Fields that convert from `i32` to other integer types now check whether
    they are in range.

  * writer: Avoid casts that may truncate.

    Fields that convert to `i32` from other integer types now check whether
    they are in range.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * Initial release.
