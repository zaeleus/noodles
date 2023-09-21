# Changelog

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
