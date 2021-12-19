# Changelog

## Unreleased

### Added

  * cram/record: Mapping quality is now stored as an `Option`.

    Valid mapping qualities are between 0 and 254, inclusive (`Some`). A
    mapping quality of 255 is considered to be missing (`None`).

## 0.9.0 - 2021-12-16

### Changed

  * cram: Update to [md-5 0.10.0].

[md-5 0.10.0]: https://crates.io/crates/md-5/0.10.0

## 0.8.3 - 2021-12-09

### Fixed

  * cram: Sync dependencies.

## 0.8.2 - 2021-12-02

### Fixed

  * cram: Require tokio's `fs` feature as a dependency ([#62]).

[#62]: https://github.com/zaeleus/noodles/issues/62

## 0.8.1 - 2021-11-18

### Fixed

  * cram: Sync dependencies.

## 0.8.0 - 2021-11-11

### Added

  * cram/crai: Add convenience write function (`crai::write`).

  * cram/crai/async: Add async writer (`crai::AsyncWriter`).

  * cram/crai/async: Add convenience write function (`crai::r#async::write`).

### Changed

  * cram: Update to Rust 2021.

## 0.7.0 - 2021-10-16

### Added

  * cram/data_container/compression_header/data_series_encoding_map/
    data_series: Add legacy TC and TN data series.

    These are no longer used in CRAM 3.0 but still need to be handled. See
    samtools/hts-specs@9a0513783826516fb8086ecf82d13631a2292f75.

  * cram/record/resolve: Handle reference skip feature in sequence resolver.

## 0.6.1 - 2021-10-02

### Fixed

  * cram: Sync dependencies.

## 0.6.0 - 2021-10-01

### Added

  * cram/reader: Add common methods to access the underlying reader: `get_ref`,
    `get_mut`, and `into_inner`.

### Fixed

  * cram/rans/decode/order_1: Fix overflow when reading frequencies.

## 0.5.1 - 2021-09-23

### Fixed

  * cram/async/reader/container/header: Fix reading starting position on the
    reference.

  * cram/async/reader/data_container/slice/header: Fix reading alignment start.

## 0.5.0 - 2021-09-19

### Added

  * cram/crai/record: Implement `Display`.

  * cram/reader: Add data container reader.

    This can be used to manually read records from slices.

  * cram/record: Add conversion to SAM record
    (`cram::Record::try_into_sam_record`).

### Changed

  * cram/record: Change alignment start to a `sam::record::Position`.

  * cram/record: Change next mate alignment start to a `sam::record::Position`.

  * cram/record/resolve: Pass compression header rather than substitution
    matrix.

    The compression header includes the substitution matrix in the preservation
    map.

### Fixed

  * cram/async/reader/data_container/slice/header: Read remainder of stream as
    optional tags.

  * cram/reader/container: Avoid casts that may truncate.

  * cram/reader/data_container/compression_header/encoding: Avoid casts that
    may truncate.

    Buffer sizes that convert from `Itf8` to `usize` now check whether they are
    in range.

  * cram/reader/data_container/slice/header: Read remainder of stream as
    optional tags.

  * cram/record/resolve: Increment feature position with operations that
    consume the read.

  * cram/record/resolve: Include last feature position.

## 0.4.0 - 2021-09-01

### Added

  * cram/async/reader: Add data container reader.

  * cram/reader: Add data container reader.

    This can be used to manually read records from slices.

### Changed

  * cram/record: `Record::read_length` is now stored as a `usize`.

### Fixed

  * cram/reader/data_container/compression_header: Avoid casts that may
    truncate.

    Buffer sizes that convert from `Itf8` to `usize` now check whether they are
    in range.

## 0.3.0 - 2021-08-19

### Added

  * cram/async: Add async header reader (`cram::AsyncReader`).

    This is a partial async CRAM reader that can only read the file definition
    and file header.

  * cram/crai/async: Add async reader (`crai::AsyncReader`).

  * cram/crai/async: Add async writer (`crai::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.2 - 2021-08-11

### Fixed

  * cram: Sync dependencies.

## 0.2.1 - 2021-07-30

### Fixed

  * cram: Sync dependencies.

## 0.2.0 - 2021-07-21

### Added

  * cram/record/tag: Add conversion from `Tag` to `sam::record::data::Field`.

### Fixed

  * cram: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * cram: Initial release.
