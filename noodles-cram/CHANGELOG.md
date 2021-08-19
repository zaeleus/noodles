# Changelog

## Unreleased

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
