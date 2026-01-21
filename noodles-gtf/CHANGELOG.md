# Changelog

## Unreleased

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.87.0.

  * gtf/record: Increase the visibility of the `attributes` module ([#360]).

    This allows access to `gtf::record::attributes::field::Value`.

[#360]: https://github.com/zaeleus/noodles/issues/360

## 0.49.0 - 2025-12-11

### Changed

  * gtf: Sync dependencies.

## 0.48.0 - 2025-11-13

### Changed

  * gtf: Sync dependencies.

## 0.47.0 - 2025-08-25

### Changed

  * gtf: Sync dependencies.

## 0.46.0 - 2025-07-12

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.45.0 - 2025-05-29

### Removed

  * gtf: Remove deprecated items.

    The following items are removed:

      * `Reader` (deprecated since 0.32.0; use `io::Reader` instead),
      * `Record::frame` (0.41.0; `Record::phase`), and
      * `Writer` (0.32.0; `io::Writer`).

## 0.44.0 - 2025-05-16

### Changed

  * gtf: Sync dependencies.

## 0.43.0 - 2025-04-13

### Changed

  * gtf: Sync dependencies.

## 0.42.0 - 2025-04-06

### Added

  * gtf/record: Implement debug formatter (`fmt::Debug`).

### Changed

  * gtf/record: Change underlying buffer to a byte string.

  * gtf/record/attributes: Unescape backslashes (`\`) and quotation marks
    (`"`).

  * gtf/io/writer/line/record/attributes/field/value: Escape backslashes (`\`)
    and quotation marks (`"`).

## 0.41.0 - 2025-03-20

### Added

  * gtf: Add lazy record (`Record`).

  * gtf/io/writer: Add feature record writer (`Writer::write_feature_record`).

### Changed

  * gtf: Replace `Record` with `gff::feature::RecordBuf`.

  * gtf: Rename `Line` to `LineBuf`.

  * gtf/io/reader: Rename line to line buf and record to record buf.

    This changes `Reader::read_line` to `Reader::read_line_buf`,
    `Reader::read_record` to `Reader::read_record_buf`, and `Reader::records`
    to `Reader::record_bufs`.

  * gtf/io/reader: Read line into `Line`.

### Deprecated

  * gtf/record: Deprecate frame usage.

    Use phase instead.

### Removed

  * gtf/line_buf: Remove parser and formatter.

    Use a deserializer (e.g., `gtf::io::Reader`) and serializer (e.g.,
    `gtf::io::Writer`), respectively, instead.

## 0.40.0 - 2025-03-08

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.39.0 - 2025-02-06

### Changed

  * gtf: Sync dependencies.

## 0.38.0 - 2025-01-24

### Fixed

  * gtf/record/attributes: Implement `AsRef<[Entry]>` ([#321]).

    This was missed in noodles-gtf 0.37.0.

[#321]: https://github.com/zaeleus/noodles/issues/321

## 0.37.0 - 2025-01-23

### Added

  * gtf/record/attributes: Add lookup by key (`Attributes::get`) ([#316]).

[#316]: https://github.com/zaeleus/noodles/issues/316

### Removed

  * gtf/record/attributes: Remove `Deref<Target = [Entry]>`.

    Use the `AsRef<[Entry]>` implementation instead.

## 0.36.0 - 2025-01-19

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.73.0.

## 0.35.0 - 2024-12-12

### Changed

  * gtf: Sync dependencies.

## 0.34.0 - 2024-11-07

### Changed

  * gtf: Sync dependencies.

## 0.33.0 - 2024-09-26

### Changed

  * gtf: Sync dependencies.

## 0.32.0 - 2024-09-09

### Added

  * gtf/io/reader: Add common methods to access the underlying I/O:
    `Reader::get_ref`, `Reader::get_mut`, and `Reader::into_inner`.

### Changed

  * gtf: Move reader (`Reader`) and writer (`Writer`) to `io` module.

### Deprecated

  * gtf: Deprecate `gtf::Reader` and `gtf::Writer`.

    Use `gtf::io::Reader` and `gtf::io::Writer`, respectively, instead.

## 0.31.0 - 2024-09-04

### Changed

  * gtf/record: Ignore trailing whitespace when parsing ([#291]).

[#291]: https://github.com/zaeleus/noodles/issues/291

### Fixed

  * gtf/record/attributes: Read string values as text between double quotes
    ([#299]).

    This allows the entry delimiter (`;`) to be used in string values.

[#299]: https://github.com/zaeleus/noodles/issues/299

## 0.30.0 - 2024-07-14

### Changed

  * gtf: Sync dependencies.

## 0.29.0 - 2024-06-17

### Changed

  * gtf: Sync dependencies.

## 0.28.0 - 2024-05-16

### Changed

  * gtf: Sync dependencies.

## 0.27.0 - 2024-05-08

### Changed

  * gtf: Sync dependencies.

## 0.26.0 - 2024-05-02

### Changed

  * gtf: Sync dependencies.

## 0.25.0 - 2024-03-28

### Changed

  * gtf: Sync dependencies.

## 0.24.0 - 2024-03-12

### Changed

  * gtf: Sync dependencies.

## 0.23.0 - 2024-01-25

### Changed

  * gtf: Sync dependencies.

## 0.22.0 - 2023-12-14

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * gtf/reader: Accept `csi::BinningIndex` for querying.

## 0.21.0 - 2023-11-14

### Changed

  * gtf: Sync dependencies.

## 0.20.0 - 2023-11-13

### Changed

  * gtf: Sync dependencies.

## 0.19.0 - 2023-10-26

### Changed

  * gtf: Sync dependencies.

## 0.18.0 - 2023-10-19

### Changed

  * gtf: Sync dependencies.

## 0.17.0 - 2023-10-12

### Changed

  * gtf: Sync dependencies.

## 0.16.0 - 2023-08-31

### Changed

  * gtf: Sync dependencies.

## 0.15.0 - 2023-08-17

### Changed

  * gtf: Sync dependencies.

## 0.14.0 - 2023-07-06

### Changed

  * gtf: Sync dependencies.

## 0.13.0 - 2023-06-29

### Changed

  * gtf: Sync dependencies.

## 0.12.0 - 2023-06-15

### Changed

  * gtf: Sync dependencies.

## 0.11.0 - 2023-06-01

### Changed

  * gtf: Sync dependencies.

## 0.10.0 - 2023-05-18

### Added

  * gtf/reader: Add query iterator (`Reader::query`) ([#158]).

[#158]: https://github.com/zaeleus/noodles/issues/158

## 0.9.0 - 2023-05-11

### Changed

  * gtf/record/attributes: Relax parser by allowing whitespace between entries
    and no final entry terminator ([#169]).

    The entry terminator (`;`) is now considered a separator rather than a
    terminator.

  * gtf/record/attributes/entry: The text representation of an entry no
    longer has a terminator (`;`).

[#169]: https://github.com/zaeleus/noodles/issues/169

## 0.8.0 - 2023-03-03

### Changed

  * gtf: Sync dependencies.

## 0.7.0 - 2023-02-03

### Added

  * gtf: Implement `std::error::Error::source` for errors.

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.64.0.

## 0.6.0 - 2022-11-18

### Added

  * gtf/writer: Add line writer (`Writer::write_line`).

## 0.5.0 - 2022-10-20

### Changed

  * gtf: Sync dependencies.

## 0.4.0 - 2022-08-16

### Changed

  * gtf: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.3.1 - 2022-06-08

### Fixed

  * gtf: Sync dependencies.

## 0.3.0 - 2022-03-29

### Changed

  * gff/record: Change start and end positions to `Position`.

## 0.2.0 - 2022-02-17

### Added

  * gtf: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.1.0 - 2021-11-11

  * gtf: Initial release.
