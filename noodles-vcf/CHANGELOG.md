# Changelog

## 0.9.0 - 2021-11-11

### Added

  * vcf/record: Add mutable getters for chromosome (`Record::chromosome_mut`),
    filters (`Record::filters_mut`), IDs (`Record::ids_mut`), info
    (`Record::info_mut`), position (`Record::position_mut`), and quality score
    (`Record::quality_score_mut`).

  * vcf/record/ids: Add ID wrapper.

  * vcf/record/info/field: Add mutable getter for value (`Field::value_mut`).

  * vcf/writer: Add common methods to access the underlying writer: `get_mut`
    and `into_inner`.

### Changed

  * vcf: Update to Rust 2021.

  * vcf/header/record/parser: Remove alphanumeric constraint on field keys.

  * vcf/record/info: `Info` no longer implements `Deref`.

    Use `AsRef<IndexMap<field::Key, Field>>` to access the lower level ordered
    map.

  * vcf/record/quality_score: Wrap only non-missing values.

    The missing quality score state is moved to `vcf::Record` as an `Option`.

### Removed

  * vcf/record/filters: Remove missing variant.

    `Filters::Missing` is removed in favor of using an `Option` in
    `vcf::Record`.

## 0.8.0 - 2021-10-16

### Added

  * vcf/record: Wrap `Genotype` field (`vcf::record::Genotypes`) ([#42]).

    This creates a new type for `Vec<vcf::record::genotypes::Genotype>`.

  * vcf/record/genotypes: Add convenience method to return a list of parsed
    genotype (`GT`) fields (`Genotypes::genotypes`) ([#42]).

  * vcf/record/genotypes/genotype: Add convenience method to parse the genotype
    (`GT`) field (`Genotype::genotype`) ([#42]).

  * vcf/record/genotypes/genotype/field/value/genotype: Add conversion from
    `Vec<Allele>` to `Genotype` ([#43]).

  * vcf/record/genotypes/genotype/field/value/genotype/allele: Add accessors
    for `position` and `phasing` ([#43]).

[#42]: https://github.com/zaeleus/noodles/issues/42
[#43]: https://github.com/zaeleus/noodles/pull/43

### Changed

  * vcf/record: Move `genotype` under `genotypes` module.

### Deprecated

  * vcf/record: Deprecated `genotype` and `Genotype` public exports.

    Use `noodles_vcf::record::genotypes::{genotype, Genotype}` instead.

## 0.7.0 - 2021-10-01

### Added

  * vcf: Increase visibility of `reader` module ([#37]).

    This allows public access to the reader iterators `Records` and `Query`.

[#37]: https://github.com/zaeleus/noodles/pull/37

## 0.6.2 - 2021-09-23

### Fixed

  * vcf: Sync dependencies.

## 0.6.1 - 2021-09-19

### Fixed

  * vcf: Sync dependencies.

## 0.6.0 - 2021-09-01

### Changed

  * vcf: Update to nom 7.0.0.

## 0.5.0 - 2021-08-19

### Changed

  * vcf: Update to tokio 1.10.0.

### Fixed

  * vcf: Define features to enable for Docs.rs.

## 0.4.0 - 2021-08-11

### Added

  * vcf/async: Add async reader (`vcf::AsyncReader`).

  * vcf/async: Add async writer (`vcf::AsyncWriter`).

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

  * vcf: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * vcf: Initial release.
