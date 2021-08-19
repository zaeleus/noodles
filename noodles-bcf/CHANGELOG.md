# Changelog

## 0.5.0 - 2021-08-19

### Changed

  * Update to tokio 1.10.0.

  * async: I/O builders are now owned/consuming builders.

    This fixes the terminal method not being able to move out of a mutable
    reference.

### Fixed

  * Define features to enable for Docs.rs.

## 0.4.0 - 2021-08-11

### Added

  * async: Add async reader (`bcf::AsyncReader`).

    This can be enabled with the `async` feature.

## 0.3.1 - 2021-08-04

### Fixed

  * Sync dependencies.

## 0.3.0 - 2021-07-30

### Added

  * header/string_map: Implement conversion from `vcf::Header`.

## 0.2.0 - 2021-07-21

### Added

  * reader: Accept any `BinningIndex` to query.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
