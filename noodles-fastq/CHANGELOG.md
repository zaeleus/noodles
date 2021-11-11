# Changelog

## 0.3.0 - 2021-11-11

### Added

  * fastq/reader: Add common methods to access the underlying reader:
    `get_ref`, `get_mut`, and `into_inner`.

### Changed

  * fastq: Update to Rust 2021.

## 0.2.0 - 2021-10-01

### Added

  * fastq/async: Add async reader (`fastq::AsyncReader`).

  * fastq/async: Add async writer (`fastq::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

### Deprecated

  * fastq/fai/record: `Record::read_name` is now `Record::name`.

  * fastq/record: `Record::read_name` is now `Record::name`.

    FASTQ record names are not necessarily read names.

## 0.1.1 - 2021-07-21

### Fixed

  * fastq: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * fastq: Initial release.
