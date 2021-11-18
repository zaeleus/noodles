# Changelog

## 0.8.0 - 2021-11-18

### Added

  * bcf/record: Implement `Debug` for `Record`.

  * bcf/record: Add getter for filters (`Record::filters`), IDs
    (`Record::ids`), genotypes (`Record::genotypes`), info (`Record::info`),
    and quality score (`Record::quality_score`).

  * bcf/record: Add mutable getters for IDs (`Record::ids_mut`), position
    (`Record::position_mut`) and quality score (`Record::quality_score_mut`).

### Changed

  * bcf/record: `bcf::Record` is no longer backed by a contiguous buffer.

    Fields are read individually when reading the record. `bcf::Record` no
    longer implements `Deref<Target = [u8]>`. `Filters`, `Info`, `Genotypes`
    now own their data.

## 0.7.0 - 2021-11-11

### Changed

  * bcf: Update to Rust 2021.

## 0.6.1 - 2021-10-16

### Fixed

  * bcf: Sync dependencies.

## 0.6.0 - 2021-10-01

### Added

  * bcf: Increase visibility of `reader` module ([#37]).

    This allows public access to the reader iterators `Records` and `Query`.

[#37]: https://github.com/zaeleus/noodles/pull/37

## 0.5.2 - 2021-09-19

### Fixed

  * bcf: Sync dependencies.

## 0.5.1 - 2021-09-01

### Fixed

  * bcf: Sync dependencies.

## 0.5.0 - 2021-08-19

### Changed

  * bcf: Update to tokio 1.10.0.

  * bcf/async: I/O builders are now owned/consuming builders.

    This fixes the terminal method not being able to move out of a mutable
    reference.

### Fixed

  * bcf: Define features to enable for Docs.rs.

## 0.4.0 - 2021-08-11

### Added

  * bcf/async: Add async reader (`bcf::AsyncReader`).

    This can be enabled with the `async` feature.

## 0.3.1 - 2021-08-04

### Fixed

  * bcf: Sync dependencies.

## 0.3.0 - 2021-07-30

### Added

  * bcf/header/string_map: Implement conversion from `vcf::Header`.

## 0.2.0 - 2021-07-21

### Added

  * bcf/reader: Accept any `BinningIndex` to query.

### Fixed

  * bcf: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * bcf: Initial release.
