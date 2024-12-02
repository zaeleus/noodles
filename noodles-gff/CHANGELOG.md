# Changelog

## Unreleased

### Added

  * gff/io/reader: Add lines iterator (`Lines`).

### Changed

  * gff: Rename `Record` to `RecordBuf.`

    This also changes `Line` to `LineBuf` and `Directive` to `DirectiveBuf`.

  * gff: Rename `lazy::Record` to `Record`.

    This also changes `lazy::Line` to `Line` and `lazy::Directive` to
    `Directive`.

  * gff/directive_buf: Change structure to hold key-optional value pairs.

    The key is now a `String`; and value, an `Option<Value>`.

  * gff/directive_buf/key: Rename name to key.

  * gff/directive_buf/key: Increase the visibilities of constants.

  * gff/directive_buf/key: Rename `START_OF_FASTA` to `FASTA`.

  * gff/directive_buf/value: Rename `Value::Other` to `Value::String`.

  * gff/io/reader: Rename line to line buf and lazy line to line.

    This changes `Lines` to `LineBufs` and `Reader::read_line` to
    `Reader::read_line_buf`, `Reader::read_lazy_line` to `Reader::read_line`.

  * gff/io/reader: Rename record to record buf.

    This changes `Reader::records` to `Reader::record_bufs`.

  * gff/io/reader/line: Skip blank lines.

    See _Generic Feature Format Version 3 (GFF3)_ (2020-08-18): "Blank lines
    should be ignored by parsers..."

  * gff/line: Hoist buffer to line.

    This moves the owned line buffer to `Line`. The structure is now a struct
    with a kind (`Kind`) field instead of an enum. The record wrapper
    (`lazy::Record`) now borrows from `Line`.

    In practice, instead of matching the line variant (e.g.,
    `Line::Record(_)`), match on the `Kind` or attempt a conversion (e.g.,
    `line.as_record()`).

  * gff/line: Wrap a directive line as `Directive<'_>`.

    This splits the line into its key-optional value components.

  * gff/record: Parse score (`Record::score`) and phase
    (`Record::phase`).

  * gff/record/attributes/field: Percent-decode components.

    Tags and values now return as `Cow<'_, str>`.

### Removed

  * gff: Remove `lazy` module.

    Types are moved to the top level, e.g., `gff::lazy::Record` is now
    `gff::Record`.

  * gff/directive_buf: Remove simple directives.

    This removes `DirectiveBuf::FeatureOntology`,
    `DirectiveBuf::AttributeOntology`, `DirectiveBuf::SourceOntology`,
    `DirectiveBuf::Species`, `DirectiveBuf::ForwardReferencesAreResolved`, and
    `DirectiveBuf::StartOfFasta`. Use `DirectiveBuf::Other` instead.

  * gff/io/reader: Remove `Reader::read_line_buf`.

    Read into a `Line` using `Reader::read_line` instead.

  * gff/line_buf: Remove formatter (`fmt::Display`).

    Use a serializer instead, e.g., `gff::io::Writer`.

### Deprecated

  * gff: Deprecate `AsyncReader`.

    Use `gff::r#async::io::Reader` instead.

## 0.39.0 - 2024-11-07

### Changed

  * gff: Sync dependencies.

## 0.38.0 - 2024-09-26

### Changed

  * gff: Sync dependencies.

## 0.37.0 - 2024-09-09

### Added

  * gff/io/reader: Add missing common methods to access the underlying I/O:
    `Reader::get_mut`, and `Reader::into_inner`.

## 0.36.0 - 2024-09-04

### Changed

  * gff: Sync dependencies.

## 0.35.0 - 2024-07-14

### Changed

  * gff: Sync dependencies.

## 0.34.0 - 2024-06-17

### Changed

  * gff: Sync dependencies.

## 0.33.0 - 2024-05-19

### Changed

  * gff: Move reader (`Reader`) and writer (`Writer`) to `io` module.

  * gff/lazy: Increase the visibility of `record` module.

  * gff/lazy/record: Parse positions (`Record::start` and `Record::end`) and
    strand (`Record::strand`).

  * gff/lazy/record/attributes: Tag the value type (`Value`).

    The raw value can still be accessed via `AsRef<str>`.

### Deprecated

  * gff: Deprecate `gff::Reader` and `gff::Writer`.

    Use `gff::io::Reader` and `gff::io::Writer`, respectively, instead.

### Removed

  * gff/lazy/record: Remove `Position` and `Strand` wrappers.

    These fields are now parsed as `noodles_core::Position` and
    `gff::record::Strand`, respectively.

## 0.32.0 - 2024-05-16

### Added

  * gff/async: Add an async reader (`gff::r#async::Reader`) ([#262]).

  * gff/lazy/line: Implement `Clone`, `Debug`, `Eq`, and `PartialEq`.

  * gff/reader: Add getter for a mutable reference to the underlying reader
    (`Reader::get_mut`).

[#262]: https://github.com/zaeleus/noodles/issues/262

### Fixed

  * gff/reader: Fix lazy reading consecutive comment or directive lines.

  * gff/reader/lazy_line: Disallow newlines to appear in fields.

    This previously allowed incomplete records to span over more than one line.

## 0.31.0 - 2024-05-08

### Changed

  * gff: Sync dependencies.

## 0.30.0 - 2024-05-02

### Changed

  * gff: Sync dependencies.

## 0.29.0 - 2024-03-28

### Changed

  * gff: Sync dependencies.

## 0.28.0 - 2024-03-12

### Changed

  * gff: Sync dependencies.

## 0.27.0 - 2024-01-25

### Changed

  * gff: Sync dependencies.

## 0.26.0 - 2023-12-14

### Changed

  * gff: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * gff/reader: Accept `csi::BinningIndex` for querying.

## 0.25.0 - 2023-11-14

### Changed

  * gff: Sync dependencies.

## 0.24.0 - 2023-11-13

### Changed

  * gff: Sync dependencies.

## 0.23.0 - 2023-10-26

### Changed

  * gff: Sync dependencies.

## 0.22.0 - 2023-10-19

### Added

  * gff/lazy: Add lazy line (`lazy::Line`) and record (`lazy::Record`).

### Changed

  * gff/record: Increase visibility of the `strand` module.

## 0.21.0 - 2023-10-12

### Changed

  * gff: Sync dependencies.

## 0.20.0 - 2023-08-31

### Changed

  * gff: Sync dependencies.

## 0.19.1 - 2023-08-24

### Fixed

  * gff/directive/name/other: Fix comparison with `&str`.

## 0.19.0 - 2023-08-17

### Changed

  * gff: Sync dependencies.

## 0.18.0 - 2023-08-03

### Added

  * gff/directive: Add support for nonstandard directives ([#190]).

    Nonstandard directives are now parsed as `Directive::Other(name::Other,
    Option<String>)`.

[#190]: https://github.com/zaeleus/noodles/issues/190

## 0.17.0 - 2023-07-20

### Added

  * gff/record/attributes/field/tag: Add standard tag constants.

### Changed

  * gff/record/attributes: Change underlying structure to an `IndexMap<Tag,
    Value>` ([#183]).

    This allows record attribute fields to have multiple values. Replace usages
    of, e.g., `attributes.iter().find(|entry| entry.key() == tag).map(|entry|
    entry.value())` with `attributes.get(tag).and_then(|values|
    values.as_string())`.

    Values are now wrapped with a typed `Value`, i.e., it can either be a
    string or array of strings. Replace `String` values with
    `Value::from(value)`.

[#183]: https://github.com/zaeleus/noodles/issues/183

## 0.16.0 - 2023-07-06

### Changed

  * gff: Sync dependencies.

## 0.15.0 - 2023-06-29

### Changed

  * gff: Sync dependencies.

## 0.14.0 - 2023-06-15

### Changed

  * gff: Sync dependencies.

## 0.13.0 - 2023-06-01

### Changed

  * gff: Sync dependencies.

## 0.12.0 - 2023-05-18

### Added

  * gff/reader: Add query iterator (`Reader::query`) ([#158]).

[#158]: https://github.com/zaeleus/noodles/issues/158

## 0.11.0 - 2023-03-03

### Changed

  * gff: Sync dependencies.

## 0.10.0 - 2023-02-03

### Added

  * gff: Implement `std::error::Error::source` for errors.

### Changed

  * gff: Raise minimum supported Rust version (MSRV) to 1.64.0.

## 0.9.0 - 2022-10-28

### Added

  * gff/writer: Add line writer (`Writer::write_line`).

## 0.8.0 - 2022-10-20

### Changed

  * gff: Sync dependencies.

## 0.7.0 - 2022-08-16

### Changed

  * gff: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.6.1 - 2022-06-08

### Fixed

  * gff: Sync dependencies.

## 0.6.0 - 2022-03-29

### Changed

  * gff/record: Change start and end positions to `Position`.

## 0.5.0 - 2022-02-17

### Added

  * gff: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.4.0 - 2021-11-11

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
