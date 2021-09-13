# Changelog

## Unreleased

### Added

  * crai/record: Implement `Display`.

  * reader: Add data container reader.

    This can be used to manually read records from slices.

  * record: Add conversion to SAM record (`cram::Record::try_into_sam_record`).

### Changed

  * record/resolve: Pass compression header rather than substitution matrix.

    The compression header includes the substitution matrix in the preservation
    map.

### Fixed

  * reader/container: Avoid casts that may truncate.

  * reader/data_container/compression_header/encoding: Avoid casts that may
    truncate.

    Buffer sizes that convert from `Itf8` to `usize` now check whether they are
    in range.

  * reader/data_container/slice/header: Read remainder of stream as optional
    tags.

  * record/resolve: Increment feature position with operations that consume the
    read.

  * record/resolve: Include last feature position.

## 0.4.0 - 2021-09-01

### Added

  * async/reader: Add data container reader.

  * reader: Add data container reader.

    This can be used to manually read records from slices.

### Changed

  * record: `Record::read_length` is now stored as a `usize`.

### Fixed

  * reader/data_container/compression_header: Avoid casts that may truncate.

    Buffer sizes that convert from `Itf8` to `usize` now check whether they are
    in range.

## 0.3.0 - 2021-08-19

### Added

  * async: Add async header reader (`cram::AsyncReader`).

    This is a partial async CRAM reader that can only read the file definition
    and file header.

  * crai/async: Add async reader (`crai::AsyncReader`).

  * crai/async: Add async writer (`crai::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.2 - 2021-08-11

### Fixed

  * Sync dependencies.

## 0.2.1 - 2021-07-30

### Fixed

  * Sync dependencies.

## 0.2.0 - 2021-07-21

### Added

  * record/tag: Add conversion from `Tag` to `sam::record::data::Field`.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
