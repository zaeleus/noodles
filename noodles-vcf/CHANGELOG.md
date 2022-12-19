# Changelog

## Unreleased

### Added

  * vcf: Implement `std::error::Error::source` for errors.

## 0.23.0 - 2022-11-29

### Added

  * vcf/header/record/key: Implement `Borrow<str>` for `Other`.

## 0.22.0 - 2022-11-18

### Added

  * vcf/header/parser: Add map value parse errors.

### Changed

  * vcf/header/record/value/map: Add actual and expected values for
    `ParseError::NumberMismatch` and `ParseError::TypeMismatch`.

### Fixed

  * vcf/header/parser: Parse record using file format ([#128]).

[#128]: https://github.com/zaeleus/noodles/issues/128

## 0.21.0 - 2022-10-28

### Changed

  * vcf: Sync dependencies.

## 0.20.0 - 2022-10-20

### Changed

  * vcf: Sync dependencies.

## 0.19.0 - 2022-09-29

### Added

  * vcf/header/record/value/map/tag: Implement `AsRef<str>` for `Tag`.

## 0.18.0 - 2022-08-16

### Added

  * vcf/header/record/value/map: Add builder (`Builder<I>`).

  * vcf/header/record/value/map/contig: Add mutable getter for length
    (`Contig::length_mut`) ([#99]).

  * vcf/header/contig/value/map/contig: Add name wrapper (`Name`).

[#99]: https://github.com/zaeleus/noodles/issues/99

### Changed

  * vcf: Raise minimum supported Rust version (MSRV) to 1.57.0.

  * vcf/header: VCF header records parsed from structured lines are now map
    values (`noodles_vcf::header::record::value::Map`).

    This moves record types (`AlternativeAllele`, `Contig`, `Filter`, `Format`,
    `Info`, and `Meta`) to map values `Map<I>`, where `I` is a specialized map
    type. A map is required to have an `ID` field; an inner `I` type that can
    have required standard fields; and can include optional fields, where the
    key is a nonstandard tag (`tag::Other`).

  * vcf/header: Change the key of `Contigs` to `contig::Name`.

  * vcf/header/record: Record holds a parsed record variant instead a key-value
    pair.

  * vcf/header/record/key: Split standard (`Standard`) and nonstandard
    (`Other`) keys.

  * vcf/header/record/value/map: Structured lines can no longer have duplicate
    field tags.

  * vcf/header/contig/value/map/contig: Change length to `usize`.

  * vcf/header/record/value/map/contig: Rename `Contig::len` to
    `Contig::length`.

### Deprecated

  * vcf/header: Deprecate `Header::records`, `Header::records_mut`,
    `header::Records`.

    Use `Header::other_records`, `Header::other_records_mut`,
    `header::OtherRecords`, respectively, instead.

### Fixed

  * vcf/header/record/value: Write surrounding angle brackets for structs.

### Removed

  * vcf/header: Remove `Pedigree`, and `Sample`.

    Use `Map<Other>` instead.

## 0.17.0 - 2022-07-05

### Changed

  * vcf: Sync dependencies.

## 0.16.1 - 2022-06-09

### Fixed

  * vcf: Sync dependencies.

## 0.16.0 - 2022-06-08

### Changed

  * vcf/record/position: Change underlying type to `usize`.

### Removed

  * vcf/record/position: Remove conversion to `i32` (`From<Position> for i32`).

    Use conversion to `usize` (`From<Position> for usize`) instead.

  * vcf/record/position: Remove fallible conversion from `i32` (`TryFrom<i32>
    for Position`).

    Use conversion from `usize` instead (`From<usize> for Position`).

  * vcf/record/position: Remove `TryFromIntError`.

    This is no longer used.

## 0.15.0 - 2022-03-29

### Changed

  * vcf: Move INFO and FORMAT keys from record field to header record.

  * vcf/record/genotypes: Add `parse` method.

    This uses FORMAT header records for type information.

  * vcf/record/genotypes/genotype: Rename `from_str_keys` to `parse`.

    `Genotype` now uses FORMAT header records for type information.

  * vcf/record/genotypes/genotype/field: Rename `from_str_key` to
    `from_str_format`.

    `Field` now uses FORMAT header records for type information.

  * vcf/record/genotypes/genotype/field/value: Rename `from_str_key` to
    `from_str_format`.

    `Value` now uses FORMAT header records for type information.

  * vcf/record/info/field/value: Rename `from_str_key` to `from_str_info`.

    `Value` now uses INFO header records for type information.

  * vcf/record/parser: Use VCF header when parsing genotypes.

### Removed

  * vcf/header: Remove record field keys.

  * vcf/record/genotypes/genotype: Remove deprecated `from_str_format`.

    Use `Genotype::parse` instead.

## 0.14.1 - 2022-03-02

### Fixed

  * vcf: Sync dependencies.

## 0.14.0 - 2022-02-17

### Added

  * vcf: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.13.0 - 2022-01-27

### Added

  * vcf/record/info: Add `clear` method to remove all fields from the info
    map.

  * vcf/record/position: Implement `Display`.

### Changed

  * vcf/record/info/field/value: Allow missing values in arrays.

    Array values (`Value::IntegerArray`, `Value::FloatArray`,
    `Value::CharacterArray`, and `Value::StringArray`) now wrap a vector of
    optional values, rather than values, e.g., `Vec<Option<i32>>` rather than
    `Vec<i32>`.

### Deprecated

  * vcf/record/genotypes/genotype: Deprecate `Genotype::from_str_format`.

    Use `Genotype::from_str_keys` instead.

### Fixed

  * vcf/record: Handle missing INFO END field value as no field.

    When a record has an INFO END field with a missing a value, i.e., `END=.`,
    it is treated as missing rather than invalid.

## 0.12.0 - 2022-01-13

### Added

  * vcf/header: Add mutable getters for `vcf::Header` fields ([#65]).

  * vcf/header: Add immutable and mutable getters for unstructured
    `vcf::Header` fields ([#69]).

  * vcf/header: Add `Records` type alias for records (`IndexMap<String,
    Vec<Record>`).

  * vcf/header/parser: Add `ParseError::StringMapPositionMismatch(actual,
    expected)`.

    This is returned when the IDX field value of a record does not match its
    relative position in a string map. It is primarily used in BCF.

  * vcf/record: Add mutable getters for genotypes, reference bases, and
    alternate bases ([#67]).

  * vcf/record/alternate_bases: Implement `DerefMut` ([#67]).

  * vcf/record/genotypes: Add mutable getter for keys ([#67]).

  * vcf/record/genotypes: Add method to return whether there are any samples
    (`Genotypes::is_empty`).

  * vcf/record/genotypes/genotype: Implement `DerefMut` ([#67]).

  * vcf/record/genotypes/keys: Implement `DerefMut` ([#67]).

  * vcf/record/position: Derive `PartialOrd` and `Ord` ([#70]).

  * vcf/record/quality_score: Derive `PartialOrd` ([#70]).

  * vcf/record/reference_bases: Implement `DerefMut` ([#67]).

[#65]: https://github.com/zaeleus/noodles/issues/65
[#67]: https://github.com/zaeleus/noodles/pull/67
[#69]: https://github.com/zaeleus/noodles/pull/69
[#70]: https://github.com/zaeleus/noodles/pull/70

### Changed

  * vcf/header/contig: Parse `IDX` field from a raw record ([#64]).

    This was previously added to the other fields map but is now a field on
    `Contig`.

### Fixed

  * vcf/header/contig: Write the IDX field value as an integer rather than a
    string ([#64]).

[#64]: https://github.com/zaeleus/noodles/issues/64

## 0.11.1 - 2021-12-09

### Fixed

  * vcf: Sync dependencies.

## 0.11.0 - 2021-12-02

### Added

  * vcf/header/format: Add conversion from `vcf::header::Record` to
    `vcf::header::Format` with a file format constraint.

  * vcf/header/info: Add conversion from `vcf::header::Record` to
    `vcf::header::Info` with a file format constraint.

  * vcf/record/genotypes: Add parser.

  * vcf/record/genotypes/keys: Implement `Default`.

### Changed

  * vcf/async/reader: `Reader::{records, query}` now parse records with the
    header.

  * vcf/reader: `Reader::{records, query}` now parse records with the header.

  * vcf/record: Move `Format` under genotypes module as `Keys`.

  * vcf/record/genotypes: Print keys on display.

  * vcf/record/genotypes/keys: Prefer key from header.

  * vcf/record/genotypes/keys: Allow empty genotypes keys.

  * vcf/record/info/field: Prefer key from header.

  * vcf/record/info/field: Allow value to be optional.

### Deprecated

  * vcf/record: Deprecate `vcf::record::Format`.

    Use `vcf::record::genotypes::Keys` instead.

### Removed

  * vcf/record/builder: Remove `Builder::{add_genotype,set_format}`.

    Use `Builder::set_genotypes` instead.

  * vcf/record: Remove `ParseError::{InvalidFormat,InvalidGenotype}`.

    These are now wrapped by `ParseError::InvalidGenotypes`.

  * vcf/record/genotypes/keys: Remove `TryFromKeyVectorError::Empty`.

    Empty genotype keys are now allowed.

  * vcf/record/genotypes: Remove `From<Vec<Genotype>>`.

    Use `Genotypes::new` instead.

## 0.10.0 - 2021-11-18

### Added

  * vcf/reader: Add common methods to access the underlying reader: `get_ref`,
    `get_mut`, and `into_inner`.

  * vcf/record/genotypes/genotype/field: Add mutable getter for value
    (`Field::value_mut`).

### Changed

  * vcf/header/format: Disable type checking when the file format is < VCF 4.3
    (#41).

  * vcf/header/info: Disable type checking when the file format is < VCF 4.3
    (#41).

[#41]: https://github.com/zaeleus/noodles/issues/41

## 0.9.0 - 2021-11-11

### Added

  * vcf/record: Add mutable getters for chromosome (`Record::chromosome_mut`),
    filters (`Record::filters_mut`), IDs (`Record::ids_mut`), info
    (`Record::info_mut`; [#52]), position (`Record::position_mut`), and quality
    score (`Record::quality_score_mut`).

  * vcf/record/ids: Add ID wrapper.

  * vcf/record/info/field: Add mutable getter for value (`Field::value_mut`;
    [#52]).

  * vcf/writer: Add common methods to access the underlying writer: `get_mut`
    and `into_inner`.

[#52]: https://github.com/zaeleus/noodles/issues/52

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
