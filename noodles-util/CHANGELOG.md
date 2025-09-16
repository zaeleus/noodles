# Changelog

## Unreleased

### Changed

  * util: Wrap inner readers and writers ([#348]).

    `alignment::io::Reader`, `alignment::io::Writer`, `variant::io::Reader`,
    and `variant::io::Writer` no longer use trait objects for the inner I/O
    object. It is now the responsibility of the caller to convert the inner I/O
    to a generic source or sink.

[#348]: https://github.com/zaeleus/noodles/pull/348

### Added

  * util/variant/async/io/reader: Add variant record reader
    (`Reader::read_record`) ([#349]).

[#349]: https://github.com/zaeleus/noodles/pull/349

## 0.70.0 - 2025-08-25

### Added

  * util/variant/io/reader: Add constructor (`Reader::new`).

    This is similar to calling `Builder::build_from_reader` with defaults.

  * util/variant/record: Implement `Debug`.

### Changed

  * util/variant/io/reader/builder: Relax reader lifetime for
    `Builder::build_from_reader`.

## 0.69.0 - 2025-07-12

### Added

  * util/variant/io/indexed_reader: Add getter for index
    (`IndexedReader::index`) ([#341]).

[#341]: https://github.com/zaeleus/noodles/issues/341

### Changed

  * util: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.68.0 - 2025-05-29

### Changed

  * util: Sync dependencies.

## 0.67.0 - 2025-05-16

### Added

  * util/alignment: Add alignment record (`alignment::Record`).

  * util/variant: Add variant record (`variant::Record`).

  * util/variant/io/reader: Add variant record reader (`Reader::read_record`).

## 0.66.0 - 2025-04-13

### Changed

  * util: Sync dependencies.

## 0.65.0 - 2025-04-06

### Changed

  * util: Sync dependencies.

## 0.64.0 - 2025-03-08

### Changed

  * util: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.63.0 - 2025-02-20

### Changed

  * util: Sync dependencies.

## 0.62.0 - 2025-02-17

### Changed

  * util: Sync dependencies.

## 0.61.0 - 2025-02-06

### Changed

  * util: Sync dependencies.

## 0.60.0 - 2025-01-24

### Changed

  * util: Sync dependencies.

## 0.59.0 - 2025-01-23

### Changed

  * util: Sync dependencies.

## 0.58.0 - 2025-01-19

### Changed

  * util: Raise minimum supported Rust version (MSRV) to 1.73.0.

## 0.57.0 - 2024-12-20

### Changed

  * util: Sync dependencies.

## 0.56.0 - 2024-12-12

### Changed

  * util: Sync dependencies.

## 0.55.0 - 2024-11-07

### Changed

  * util: Sync dependencies.

## 0.54.0 - 2024-10-22

### Changed

  * util: Sync dependencies.

## 0.53.1 - 2024-09-26

### Changed

  * util: Sync dependencies.

### Fixed

  * util/alignment/io/indexed_reader/builder: Use custom index if set when
    building from a path ([#303]).

[#303]: https://github.com/zaeleus/noodles/issues/303

## 0.53.0 - 2024-09-12

### Changed

  * util: Sync dependencies.

## 0.52.0 - 2024-09-09

### Changed

  * util: Sync dependencies.

## 0.51.0 - 2024-09-04

### Added

  * util/alignment: Add async reader and writer
    (`alignment::r#async::io::Reader` and `alignment::r#async::io::Writer`)
    ([#286]).

[#286]: https://github.com/zaeleus/noodles/issues/286

## 0.50.0 - 2024-08-04

### Added

  * util/variant: Add async reader (`variant::r#async::io::Reader`) and writer
    (`variant::r#async::io::Writer`).

## 0.49.0 - 2024-07-14

### Changed

  * util: Sync dependencies.

## 0.48.0 - 2024-06-17

### Changed

  * util: Sync dependencies.

## 0.47.0 - 2024-06-06

### Changed

  * util: Sync dependencies.

## 0.46.0 - 2024-05-31

### Changed

  * util: Sync dependencies.

## 0.45.0 - 2024-05-16

### Changed

  * util: Sync dependencies.

## 0.44.0 - 2024-05-08

### Changed

  * util: Sync dependencies.

## 0.43.0 - 2024-05-02

### Changed

  * util: Sync dependencies.

## 0.42.0 - 2024-04-22

### Changed

  * util: Sync dependencies.

## 0.41.0 - 2024-04-11

### Changed

  * util: Sync dependencies.

## 0.40.0 - 2024-04-04

### Changed

  * util/variant: Move readers (`Reader` and `IndexedReader`) and writer
    (`Writer`) to `io` module.

  * util/variant/io/indexed_reader: Change iterators (`IndexedReader::records`
    and IndexedReader::query`) to return an iterator over
    `vcf::variant::Record` trait objects.

  * util/variant/io/reader: Change `Reader::records` to return an iterator
    over `vcf::variant::Record` trait objects.

  * util/variant/io/writer: Change `Writer::write_record` to accept `&dyn
    vcf::variant::Record`.

### Fixed

  * util/alignment/io/indexed_reader: Use format record iterators in records
    iterator (`IndexedReader::records`).

## 0.39.0 - 2024-03-28

### Changed

  * util: Sync dependencies.

## 0.38.0 - 2024-03-12

### Changed

  * util: Sync dependencies.

## 0.37.0 - 2024-02-22

### Changed

  * util: Sync dependencies.

## 0.36.0 - 2024-02-15

### Changed

  * util: Sync dependencies.

## 0.35.0 - 2024-02-08

### Changed

  * util: Sync dependencies.

## 0.34.1 - 2024-02-04

### Fixed

  * util: Sync dependencies.

## 0.34.0 - 2024-02-01

### Changed

  * util: Sync dependencies.

## 0.33.0 - 2024-01-25

### Added

  * util/alignment: Add `iter` module from `noodles_sam::alignment`.

    This includes a pileup iterator to calculate sequence depth.

### Changed

  * util/alignment: Move readers (`Reader` and `IndexedReader`) and writer
    (`Writer`) to `io` module.

  * util/alignment/io/indexed_reader: Change `IndexedReader::query` to return
    an iterator over `io::Result<Box<dyn sam::alignment::Record>>`.

  * util/alignment/io/reader: Change `Reader::records` to return an iterator
    over `io::Result<Box<dyn sam::alignment::Record>>`.

  * util/alignment/io/writer: Change `Writer::write_record` to accept
    `&dyn sam::alignment::Record`.

### Removed

  * util/variant: Remove `Compression`.

    This was deprecated in noodles-util 0.20.0 Use `CompressionMethod` instead.

## 0.32.0 - 2023-12-14

### Changed

  * util: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.31.0 - 2023-11-14

### Changed

  * util: Sync dependencies.

## 0.30.0 - 2023-11-13

### Changed

  * util: Sync dependencies.

## 0.29.0 - 2023-11-02

### Changed

  * util: Sync dependencies.

## 0.28.0 - 2023-10-26

### Changed

  * util: Sync dependencies.

## 0.27.0 - 2023-10-23

### Changed

  * util: Sync dependencies.

## 0.26.0 - 2023-10-19

### Changed

  * util: Sync dependencies.

## 0.25.0 - 2023-10-12

### Changed

  * util: Sync dependencies.

## 0.24.0 - 2023-09-21

### Changed

  * util: Sync dependencies.

## 0.23.0 - 2023-09-14

### Changed

  * util: Sync dependencies.

## 0.22.0 - 2023-08-31

### Changed

  * util: Sync dependencies.

## 0.21.0 - 2023-08-24

### Changed

  * util: Sync dependencies.

## 0.20.0 - 2023-08-17

### Added

  * util/alignment: Add an indexed reader (`alignment::IndexedReader`).

  * util/alignment: Add `CompressionMethod` to override I/O compression.

    The compression method can now be overridden for most formats.

  * util/alignment/reader/builder: Add support for raw BAM.

### Changed

  * util/alignment/writer/builder: `Builder::build_from_writer` now returns
    an `io::Result<Writer>`.

### Deprecated

  * util/variant: Deprecate `Compression`.

    Use `CompressionMethod` instead.

  * util/variant: Deprecate `variant::*::Builder::set_compression`.

    Use `variant::*::Builder::set_compression_method` instead.

### Removed

  * util/alignment/format: Remove `Format::SamGz`.

    Set the format-compression method to `(Format::Sam,
    Some(CompressionMethod::Bgzf))` instead.

## 0.19.1 - 2023-08-04

### Fixed

  * util/variant/indexed_reader/builder: Avoid discarding buffer after file
    type detection.

## 0.19.0 - 2023-08-03

### Added

  * util/variant/indexed_reader/builder: Add compression
    (`Builder::set_compression`), format (`Builder::set_format`), and index
    setters (`Builder::set_index`).

  * util/variant/indexed_reader/builder: Add `Builder::build_from_reader`.

## 0.18.0 - 2023-07-27

### Added

  * util/variant: Add an indexed reader (`variant::IndexedReader`) ([#188]).

[#188]: https://github.com/zaeleus/noodles/issues/188

## 0.17.0 - 2023-07-20

### Changed

  * util: Sync dependencies.

## 0.16.0 - 2023-07-06

### Changed

  * util: Sync dependencies.

### Fixed

  * util/alignment/variant/builder: Fix format detection when the first BGZF
    block is incomplete.

## 0.15.0 - 2023-06-29

### Changed

  * util: Sync dependencies.

### Fixed

  * util/alignment/reader/builder: Fix format detection when the first BGZF
    block is incomplete (#179).

[#179]: https://github.com/zaeleus/noodles/issues/179

## 0.14.0 - 2023-06-15

### Changed

  * util: Sync dependencies.

## 0.13.0 - 2023-06-08

### Changed

  * util: Sync dependencies.

## 0.12.0 - 2023-06-01

### Changed

  * util: Sync dependencies.

## 0.11.0 - 2023-05-18

### Changed

  * util: Sync dependencies.

## 0.10.0 - 2023-05-11

### Added

  * util/alignment/format: Add support for bgzipped SAM (`*.sam.gz`).

## 0.9.0 - 2023-05-04

### Changed

  * util: Sync dependencies.

## 0.8.0 - 2023-04-27

### Changed

  * util: Sync dependencies.

## 0.7.0 - 2023-04-06

### Changed

  * util: Sync dependencies.

## 0.6.0 - 2023-03-14

### Added

  * util: Add variant writer (`variant::Writer`) ([#150]).

    This is a high-level writer that abstracts writing both VCF and BCF. It can
    autodetect the output format and compression type at runtime.

[#150]: https://github.com/zaeleus/noodles/issues/150

## 0.5.0 - 2023-03-03

### Added

  * util: Add variant reader (`variant::Reader`) ([#149]).

    This is a high-level reader that abstracts reading both VCF and BCF. It can
    autodetect the input format and compression type at runtime.

[#149]: https://github.com/zaeleus/noodles/pull/149

## 0.4.0 - 2023-02-03

### Changed

  * util: Raise minimum supported Rust version (MSRV) to 1.64.0.

## 0.3.1 - 2022-11-29

### Changed

  * util: Sync dependencies.

## 0.3.0 - 2022-11-18

### Changed

  * util/alignment/reader/builder: `Builder::build_from_reader` is no longer
    constrained to `Seek` ([#130]).

[#130]: https://github.com/zaeleus/noodles/issues/130

## 0.2.0 - 2022-10-28

### Changed

  * util: Sync dependencies.

## 0.1.0 - 2022-10-20

  * util: Initial release.
