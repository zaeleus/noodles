# Changelog

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
