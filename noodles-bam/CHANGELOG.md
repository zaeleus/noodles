# Changelog

## Unreleased

### Added

  * bai/index: Implemented `BinningIndex` for `Index`.

  * bai/index: Added `query` method to find chunks that intersect the given
    region.

  * bai/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

  * reader: Accept any `BinningIndex` to query.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

  * reader: Avoid casts that may truncate.

    Fields that convert to `u32` from other integer types now check whether
    they are in range.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * bai/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * Initial release.
