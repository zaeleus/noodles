# Changelog

## Unreleased

### Added

  * async: Add async reader (`sam::AsyncReader`).

    The async reader can be enabled with the `async` feature.

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
