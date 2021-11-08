# Changelog

## Unreleased

### Changed

  * gff: Update to Rust 2021.

## 0.3.0 - 2021-10-16

### Added

  * gff/record: Implement `fmt::Display`.

    The string representation of a `Record` is its serialized tabular form.

### Changed

  * gff/record: Disallow reference sequence names to start with '>'.

## 0.2.0 - 2021-09-19

### Added

  * gff/record/attributes/entry: Accept `Into<String>` for key and value.

### Changed

  * gff/record/attributes/entry: Return `ParseError::Invalid` when no `=`
    separator is present.

    This previously would return `ParseError::MissingValue`.

## 0.1.1 - 2021-07-21

### Fixed

  * gff: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * gff: Initial release.
