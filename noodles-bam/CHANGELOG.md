# Changelog

## 0.2.1 - 2021-07-30

### Fixed

  * bai/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * bai/index: Implemented `BinningIndex` for `Index`.

  * bai/index: Added `query` method to find chunks that intersect the given
    region.

  * bai/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

  * reader: Accept any `BinningIndex` to query.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

  * bai/reader: Avoid casts that may truncate.

    Fields that convert from `u32` to other integer types now check whether
    they are in range.

  * bai/writer: Avoid casts that may truncate.

    Fields that convert to `u32` from other integer types now check whether
    they are in range.

  * writer/record: Fix bin calculation start position.

    This incorrectly used a 1-based position rather than 0-based.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * bai/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * Initial release.
