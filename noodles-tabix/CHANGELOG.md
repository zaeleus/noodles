# Changelog

## Unreleased

### Added

  * index: Implemented `BinningIndex` for `Index`.

  * index/reference_sequence: Implemented `BinningIndexReferenceSequence` for
    `ReferenceSequence`.

### Deprecated

  * index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

### Fixed

  * writer: Avoid casts that may truncate.

    Fields that convert to `i32` from other integer types now check whether
    they are in range.

### Removed

  * index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * Initial release.
