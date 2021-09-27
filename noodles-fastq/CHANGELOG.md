# Changelog

## Unreleased

### Added

  * async: Add async reader (`fastq::AsyncReader`).

  * async: Add async writer (`fastq::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

### Deprecated

  * fai/record: `Record::read_name` is now `Record::name`.

  * record: `Record::read_name` is now `Record::name`.

    FASTQ record names are not necessarily read names.

## 0.1.1 - 2021-07-21

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
