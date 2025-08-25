# Changelog

## 0.28.0 - 2025-08-25

### Changed

  * bed/io/writer/record: Constrain feature positions to 64-bit unsigned
    integers.

## 0.27.0 - 2025-07-12

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.26.0 - 2025-05-29

### Removed

  * bed: Remove deprecated items.

    The following items are removed:

      * `Reader` (0.14.0; use `io::Reader` instead) and
      * `Writer` (0.14.0; `io::Writer`).

## 0.25.0 - 2025-05-16

### Changed

  * bed: Sync dependencies.

## 0.24.0 - 2025-04-13

### Changed

  * bed: Sync dependencies.

## 0.23.0 - 2025-04-06

### Changed

  * bed: Sync dependencies.

## 0.22.0 - 2025-03-08

### Added

  * bed/fs: Add file indexer.

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.21.0 - 2025-02-06

### Changed

  * bed: Sync dependencies.

## 0.20.0 - 2025-01-23

### Changed

  * bed: Sync dependencies.

## 0.19.0 - 2025-01-19

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.73.0.

## 0.18.0 - 2024-12-12

### Changed

  * bed: Sync dependencies.

## 0.17.0 - 2024-09-26

### Changed

  * bed: Update to lexical-core 1.0.0.

## 0.16.0 - 2024-09-04

### Added

  * bed: Add a record view (`Record`).

  * bed/feature: Add `Record` trait to represent an opaque feature record.

    noodles-bed now only supports BED3+, BED4+, BED5+, and BED6+. Records with
    more standard fields should use BED6+, e.g., BED12 can be represented as
    BED6+6.

  * bed/feature/other_fields: Values can now be typed (`Value`).

    This allows setting values other than strings.

  * bed/io: Add builders (`reader::Builder` and `writer::Builder`).

### Changed

  * bed: Move `Record` to `feature::RecordBuf`.

  * bed/record: Represent a feature end of 0 as `None`.

    This changes `Record::feature_end` from `Position` to `Option<Position>`,
    which allows modeling features before the record's reference sequence.

  * bed/feature: Move `record_buf::Strand` to `record::Strand`.

    `Strand` is now a simple representation of the forward and reverse options.
    It no longer parses or serializes the field.

  * bed/io: `Reader` and `Writer` must be annotated with the number of expected
    standard fields.

    `Reader` and `Writer` now require the number of expected standard fields as
    a const generic argument, e.g., for BED3+ records, use `Reader<3, _>` and
    `Writer<3, _>`.

    This also means that dynamically sized records are no longer supported.

### Removed

  * bed/feature/record_buf: Remove `str::FromStr` and `fmt::Display`.

    Use a BED record reader and writer, respectively, instead.

  * bed/feature/record_buf: Remove `Color`.

    This is no longer used.

  * bed/io/reader: Remove `Reader::read_line`.

    Read into a `Record<N>` buffer instead.

## 0.15.0 - 2024-06-17

### Changed

  * bed/record/name: Remove validation on parse ([#271]).

    This allows for more relaxed input in memory and is checked on
    serialization instead.

[#271]: https://github.com/zaeleus/noodles/issues/271

### Fixed

  * bed/record: Fix parse error kind returned from an invalid score ([#270]).

[#270]: https://github.com/zaeleus/noodles/pull/270

## 0.14.0 - 2024-06-06

### Changed

  * bed: Move reader (`Reader`) and writer (`Writer`) to `io` module.

### Deprecated

  * bed: Deprecate `bed::Reader` and `bed::Writer`.

    Use `bed::io::Reader` and `bed::io::Writer`, respectively, instead.

## 0.13.0 - 2024-05-08

### Changed

  * bed: Sync dependencies.

## 0.12.0 - 2024-01-25

### Changed

  * bed: Sync dependencies.

## 0.11.0 - 2023-12-14

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.10.0 - 2023-06-15

### Changed

  * bed: Sync dependencies.

## 0.9.0 - 2023-05-18

### Added

  * bed/record/score: Add a const constructor (`Score::new`), getter
    (`Score::get`), and min (`Score::MIN`) and max (`Score::MAX`) constants.

### Fixed

  * bed/record/name: Fix valid name character set.

    This previously disallowed the tilde (`~`) character.

## 0.8.0 - 2023-03-03

### Changed

  * bed: Sync dependencies.

## 0.7.0 - 2023-02-03

### Added

  * bed/record: Implement `std::error::Error::source` for errors.

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.64.0.

## 0.6.0 - 2022-11-18

### Added

  * bed/record/color: Add getters for color components (`Color::red`,
    `Color::green`, and `Color::blue`).

## 0.5.0 - 2022-10-20

### Changed

  * bed: Sync dependencies.

## 0.4.0 - 2022-08-16

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.3.0 - 2022-06-08

### Added

  * bed/record: Add support for BED7+ (`Record<7>`), BED8+ (`Record<8>`), BED9+
    (`Record<9>`), and BED12+ (`Record<12>`).

## 0.2.0 - 2022-03-29

### Changed

  * bed/record: Wrap start and end positions.

    Positions are now 1-based, inclusive.

## 0.1.0 - 2022-02-17

  * bed: Initial release.
