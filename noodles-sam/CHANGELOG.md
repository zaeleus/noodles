# Changelog

## Unreleased

### Added

  * sam/record: Add mutable getters for flags (`Record::flags_mut`), read name
    (`Record::read_name_mut`), mapping quality (`Record::mapping_quality_mut`),
    and template length (`Record::template_length_mut`).

## 0.4.0 - 2021-09-23

### Changed

  * async/reader: Handle CRLF newlines and missing final newline.

  * header/record: Require delimiter split when parsing.

    This ensures the delimiter exists when tokenizing. It removes
    `ParseError::MissingValue` for `ParseError::Invalid` and
    `ParseError::MissingTag` for `ParseError::InvalidField`.

  * header/record: Validate field values (`/[ -~]+/`).

  * reader: Handle CRLF newlines and missing final newline.

## 0.3.0 - 2021-09-19

### Added

  * async: Add async reader (`sam::AsyncReader`).

  * async: Add async writer (`sam::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

  * reader: Add method to return the underlying reader (`Reader::into_inner`).

  * writer: Add method to return the underlying writer (`Writer::into_inner`).

### Changed

  * header/program/builder: Return error from `build`.

    This previously panicked if the ID was not set.

  * header/read_group/builder: Return error from `build`.

    This previously panicked if the ID was not set.

## 0.2.2 - 2021-08-19

### Fixed

  * Sync dependencies.

## 0.2.1 - 2021-08-11

### Fixed

  * Sync dependencies.

## 0.2.0 - 2021-07-30

### Added

  * header/reference_sequence: Add alternative locus (`AH`) parser.

### Changed

  * header/reference_sequence: Parse alternative locus (`AH`) value.

    This is no longer stored as a raw `String`.

## 0.1.1 - 2021-07-21

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
