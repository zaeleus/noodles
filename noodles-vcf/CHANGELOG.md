# Changelog

## Unreleased

### Added

  * record/genotype/field/value/genotype: Add conversion from `Vec<Allele>` to
    `Genotype` ([#43]).

  * record/genotype/field/value/genotype/allele: Add accessors for `position`
    and `phasing` ([#43]).

[#43]: https://github.com/zaeleus/noodles/pull/43

## 0.7.0 - 2021-10-01

### Added

  * Increase visibility of `reader` module ([#37]).

    This allows public access to the reader iterators `Records` and `Query`.

[#37]: https://github.com/zaeleus/noodles/pull/37

## 0.6.2 - 2021-09-23

### Fixed

  * Sync dependencies.

## 0.6.1 - 2021-09-19

### Fixed

  * Sync dependencies.

## 0.6.0 - 2021-09-01

### Changed

  * Update to nom 7.0.0.

## 0.5.0 - 2021-08-19

### Changed

  * Update to tokio 1.10.0.

### Fixed

  * Define features to enable for Docs.rs.

## 0.4.0 - 2021-08-11

### Added

  * async: Add async reader (`vcf::AsyncReader`).

  * async: Add async writer (`vcf::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.3.0 - 2021-08-04

### Added

  * vcf/header/alternative_allele: Accept `Into<String>` for description in
    constructor.

## 0.2.0 - 2021-07-30

### Added

  * vcf/header/contig: Accept `Into<String>` arguments in constructor.

  * vcf/header/filter: Accept `Into<String>` arguments in constructor.

## 0.1.1 - 2021-07-21

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
