# Changelog

## Unreleased

## Deprecated

  * fasta: Deprecate async re-exports (`AsyncReader` and `AsyncWriter`).

    Use `fastq::r#async::io::Reader` and `fastq::r#async::io::Writer` instead.

## 0.15.0 - 2024-09-26

### Changed

  * fastq/fai: Move reader (`Reader`) and writer (`Writer`) to `io` module.

### Deprecated

  * fastq/fai: Deprecate `Reader` and `Writer`.

    Use `fai::io::Reader` and `fai::io::Writer`, respectively, instead.

## 0.14.0 - 2024-08-04

### Added

  * fastq/io/reader/record/definition: Support horizontal tabs (`\t`) as a
    definition separator ([#287]).

    The reader now splits the definition name and description on either a space
    (` `) or horizontal tab (`\t`).

  * fastq/io/writer: Add builder.

  * fastq/io/writer: Support setting a custom definition separator
    (`Builder::set_definition_separator`) ([#287]).

    The default definition separator is space (` `).

  * fastq/io/writer: Add additional methods to get the underlying writer
    (`Writer::get_mut` and `Writer::into_inner`).

[#287]: https://github.com/zaeleus/noodles/issues/287

### Changed

  * fastq/record/definition: Change underlying types to byte strings
    (`BString`).

## 0.13.0 - 2024-07-14

### Removed

  * fastq/record: Remove `fmt::Display`.

    Write a record using `io::Writer` instead.

## 0.12.0 - 2024-06-17

### Added

  * fastq/fai: Add common methods to access the underlying I/O.

## 0.11.0 - 2024-05-31

### Added

  * fastq/async/io/reader: Add common methods to access the underlying reader:
    `Reader::get_ref`, `Reader::get_mut`, and `Reader::into_inner`.

  * fastq/async/io/writer: Add `Writer::get_mut`.

### Changed

  * fastq: Move reader (`Reader`), writer (`Writer`), and indexer (`Indexer`)
    to `io` module.

### Deprecated

  * fastq: Deprecate `fastq::Reader`, `fastq::Writer`, `fastq::Indexer`, and
    `fasta::index`.

    These are moved under the `io` module.

## 0.10.0 - 2023-12-14

### Changed

  * fastq: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.9.0 - 2023-10-12

### Deprecated

  * fastq/fai/record: Deprecate `Record::len`.

    Use `Record::length` instead.

## 0.8.0 - 2023-05-18

### Changed

  * fastq/fai/record: Accept `Into<String>` for the record name.

  * fastq/record: Add `Definition` wrapper ([#165]).

    Like `fasta::record::Definition`, `fastq::record::Definition` wraps the
    record name and optional description. Change usages of
    `fastq::Record::new(name, sequence, quality_scores)` to
    `fastq::Record::new(Definition::new(name, description), sequence,
    quality_scores)`.

[#165]: https://github.com/zaeleus/noodles/issues/165

## 0.7.1 - 2023-05-11

### Fixed

  * fastq/reader/record/definition: Fix panic on trailing newline removal
    ([#166]).

    This previously checked the wrong buffer when removing the trailing
    newline.

[#166]: https://github.com/zaeleus/noodles/issues/166

## 0.7.0 - 2023-04-27

### Changed

  * fastq/indexer: Trim description from read name.

    The index record no longer includes the description from the read name
    line.

  * fastq/reader: The description is now considered part of the name line.

    This is similar to the FASTA definition, where the FASTQ definition is
    "@<name> <description>". The content of the "plus" line is discarded.

  * fastq/writer: Write description as part of name line ([#161]).

    This is no longer written to the "plus" line.

[#161]: https://github.com/zaeleus/noodles/issues/161

## 0.6.0 - 2023-02-03

### Added

  * fastq/fai/record: Implement `std::error::Error::source` for `ParseError`.

### Changed

  * fastq: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * fastq: Remove `Record::name`.

    This was deprecated in noodles-fastq 0.2.0. Use `Record::name`
    instead.

## 0.5.1 - 2022-10-28

### Fixed

  * fastq/reader: Increase the visibility of the module (`reader`) ([#118]).

    This allows access to the `reader::Records` iterator.

[#118]: https://github.com/zaeleus/noodles/issues/118

## 0.5.0 - 2022-02-17

### Added

  * fastq: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.4.0 - 2022-01-27

### Added

   * fastq/record: Add description field (`Record::description`).

   * fastq/record: Add mutable getters for name (`Record::name_mut`),
     sequence (`Record::sequence_mut`), description
     (`Record::description_mut`), and quality scores
     (`Record::quality_scores_mut`).

### Changed

  * fastq/async/reader: Ensure the record description (line 3) is prefixed
    with a plus sign (`+`).

  * fastq/reader: Ensure the record description (line 3) is prefixed with a
    plus sign (`+`).

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
