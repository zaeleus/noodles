# Changelog

## Unreleased

### Added

  * vcf/header/record/key/other: Implement `FromStr` for `Other`.

  * vcf/header/record/value/map/contig: Add getter (`Map::<Contig>::url`) and
    mutable getter (`Map::<Contig>::url_mut`) for URL.

  * vcf/record/info: Change getters to accept keys that implement
    `Equivalent<Key>`.

### Changed

  * vcf/header: Change other records to a collection of either unstructured or
    structured records.

    This changes `OtherRecords` to a be `IndexMap<key::Other,
    value::Collection>`, where the collection can be a list of strings or a map
    of `Map<Other>` values. Change usages of `header::record::value::Other` to
    `header::record::Value`.

  * vcf/header/format: Move `Key` under `record::genotypes::keys` module.

  * vcf/header/info: Move `Key` under `record::info::field` module.

  * vcf/header/record/value/map: Increase the visibility of `Tag` for inner
    maps.

  * vcf/header/record/value/map/builder: Accept `Into<String>` for other field
    values (`Builder::insert`).

  * vcf/header/record/value/map/other: Move record ID from record to map key.

    `Map<Other>` no longer holds the record ID. Use the collection key instead.

  * vcf/indexed_reader/builder: Attempt to load either an associated tabix
    index (`<src>.tbi`) or CSI (`<src>.csi`), in that order.

  * vcf/record/alternate_bases: Change `Deref` target from `[Allele]` to
    `Vec<Allele>`.

### Fixed

  * vcf/header/record/value/map/builder: Implement builder traits for generic
    builders.

  * vcf/header/record/value/map/contig: Include URL in `Display`
    implementation.

### Removed

  * vcf/header/format: Remove `format` module.

    For genotype keys, use `record::genotypes::keys::Key` instead.

  * vcf/header/info: Remove `info` module.

    For info keys, use `record::info::field::Key` instead.

## 0.27.0 - 2023-04-06

### Added

  * vcf/header/record/value/map/contig: Add `md5` and `URL` as a standard tags.

  * vcf/indexed_reader: Add getter for index (`IndexedReader::index`).

  * vcf/record: Implement `Default`.

    This is meant as a convenience to more easily create a mutable buffer.

  * vcf/record: Implement `TryFrom<(&Header, &str)>`.

    This can be used to parse a raw record.

  * vcf/record/ids: Implement `Extend<Id>` and `FromIterator<Id>`.

### Changed

  * vcf/async/reader: Change `Reader::read_header` to return a parsed header
    (`Header`).

    This no longer returns a raw string.

  * vcf/async/reader: Change `Reader::read_record` to read into a `Record`
    buffer.

    This now requires a `Header` and `Record` buffer. It no longer returns the
    raw string.

  * vcf/reader: Change `Reader::read_header` to return a parsed header
    (`Header`).

    This no longer returns a raw string.

  * vcf/reader: Change `Reader::read_record` to read into a `Record` buffer.

    This now requires a `Header` and `Record` buffer. It no longer returns the
    raw string.

  * vcf/record/chromosome: Allow "." as a valid name.

  * vcf/record/genotypes: Own all sample values.

    `Genotypes` now owns all sample values, rather than `Sample`. Genotype
    values changed to a list of a list of optional values
    (`Vec<Vec<Option<Value>>>`).

  * vcf/record/genotypes: Rename `Genotype` to `Sample`.

  * vcf/record/info: Move missing field value parsing to record parser.

    `Info` no longer handles "." as missing.

  * vcf/record/parser: Use `vcf::reader::parse_record` to parse records.

  * vcf/record/reference_bases: Move missing field value parsing to record
    parser.

    `ReferenceBases` no longer handles "." as missing.

  * vcf/writer: Require header when writing a record (`Writer::write_record`).

### Deprecated

  * vcf/record: Deprecate `Record::try_from_str`.

    Use `TryFrom<(&Header, &str)>` instead.

### Removed

  * vcf/record: Remove `ParseError`.

    Use `vcf::reader::record::ParseError` instead.

  * vcf/record: Remove `Field`.

    This is no longer used.

  * vcf/record/genotypes: Remove `Deref` and `DerefMut` for `Genotypes`.

  * vcf/record/genotypes: Remove `Values`.

    Use `Sample` instead.

  * vcf/record/genotypes/values: Remove `field` module.

    This removes the unused `Field` struct. The `values::field::value` module
    is moved to `sample::value`.

## 0.26.0 - 2023-03-14

### Added

  * vcf: Add an indexed reader (`vcf::IndexedReader`).

  * vcf/reader: Add builder (`vcf::reader::Builder`).

    The builder is able to construct a reader from a path
    (`Builder::build_from_path`), which can open raw VCF files (`*.vcf`)
    and bgzipped VCF (`*.vcf.gz`) files.

  * vcf: Add a variant writer trait (`VariantWriter`) ([#150]).

    This is a generalization for writing VCF-like variant formats.

  * vcf/writer: Implement `VariantWriter` ([#150]).

[#150]: https://github.com/zaeleus/noodles/issues/150

### Changed

  * vcf/record/alternate_bases: Move missing field value parsing to record
    parser.

    `AlternateBases` no longer handles "." as missing.

### Fixed

  * vcf/record/alternate_bases/allele: Ensure breakend has more than one
    character.

## 0.25.0 - 2023-03-03

### Added

  * vcf: Add a variant reader trait (`VariantReader`) ([#149]).

    This is a generalization for reading VCF-like variant formats.

  * vcf/header: Add header parser (`header::Parser`).

    This can be used to customize how to parse the header.

  * vcf/header/format/key: Implement `Borrow<str>`.

  * vcf/header/format/key: Add VCF 4.4 reserved keys.

  * vcf/header/info/key: Implement `Borrow<str>`.

  * vcf/header/info/key: Add VCF 4.4 reserved keys.

  * vcf/reader: Implement `VariantReader` ([#149]).

  * vcf/record/info/field/value: Implement `TryFrom<(Number, Type, &str)>`.

  * vcf/record/genotypes/genotype/field/value: Implement `TryFrom<(Number,
    Type, &str)>`.

[#149]: https://github.com/zaeleus/noodles/pull/149

### Changed

  * vcf/header: Move `header::format::Type` to record map value.

    Use `header::record::value::map::format::Type` instead.

  * vcf/header: Move `header::info::Type` to record map value.

    Use `header::record::value::map::info::Type` instead.

  * vcf/header/format/key: Split standard (reserved) and other (non-reserved)
    keys.

    Change usages of, e.g, `Key::Genotype` to `key::GENOTYPE`.

  * vcf/header/info/key: Split standard (reserved) and other (non-reserved)
    keys.

    Change usages of, e.g, `Key::TotalDepth` to `key::TOTAL_DEPTH`.

  * vcf/header/record/parser: Disallow empty values for unstructured lines.

  * vcf/record/genotypes/genotype/field/value/genotype: Infer phasing of first
    allele.

  * vcf/record/genotypes/genotype/field/value/genotype/allele: The phasing is
    now required.

  * vcf/record/ids: Move missing field value parsing to record parser.

    `Ids` no longer handles "." as missing.

  * vcf/record/ids/id: Disallow `.` as a valid identifier.

### Fixed

  * vcf/reader: Fix infinite loop when an input is only a header with no final
    newline.

    This makes the header line reader consistent with the record line reader.

  * vcf/record/genotypes/genotype: Fail parsing when there are more values than
    keys.

    This previously would silently drop values after the last key.

## 0.24.0 - 2023-02-03

### Added

  * vcf: Implement `std::error::Error::source` for errors.

  * vcf/record/genotypes/genotype: Implement `Extend<(Key, Option<Value>)>`.

  * vcf/record/genotypes/genotype: Implement `FromIterator<(Key,
    Option<Value>)>`.

  * vcf/record/info: Implement `Extend<(Key, Option<Value>)>`.

  * vcf/record/info: Implement `FromIterator<(Key, Option<Value>)>`.

### Changed

  * vcf: Raise minimum supported Rust version (MSRV) to 1.64.0.

  * vcf/header: Change inserting `OtherRecords` as an `Key::Other`-`Value`
    pair.

  * vcf/header: Change `Header::get` to accept keys that implement
    `Equivalent<Other>`.

  * vcf/header/builder: Change inserting ALT, contig, FILTER, FORMAT, INFO, and
    META records as an ID-map pair.

  * vcf/header/record/value/map: Move records IDs from record to map key.

    `Map<I>` no longer holds the record ID. Use the collection key instead.

  * vcf/header/record/value/map/format: Change conversion from key value to a
    reference (`&format::Key`).

  * vcf/header/record/value/map/info: Change conversion from key value to a
    reference (`&info::Key`).

  * vcf/record/info: Change underlying map to `Key`-`Option<Value>` rather than
    `Key`-`Field`.

  * vcf/record/genotypes/genotype: Change underlying map to
    `Key`-`Option<Value>` rather than `Key`-`Field`.

  * vcf/record/genotypes/genotype: Change `TryFrom<Vec<Field>>` to
    `TryFrom<Vec<(Key, Option<Value>)>>`.

### Removed

  * vcf/header: Remove `Records`.

    This was deprecated in noodles-vcf 0.18.0. Use `OtherRecords`
    instead.

  * vcf/header: Remove `Header::records` and `Header::records_mut`.

    These were deprecated in noodles-vcf 0.18.0. Use
    `Header::other_records` and `Header::other_records_mut`,
    respectively, instead.

  * vcf/record: Remove re-exports of `record::genotypes::genotype`,
    `record::genotypes::Genotype`, and `record::genotypes::Keys as
    Format`.

    These were deprecated in 0.11.0.

  * vcf/record/info: Remove `TryFrom<Vec<Field>>` for `Info`.

    Insert fields into the map using `Info::insert` instead.

  * vcf/record/info/field: Remove `Field`.

    Use `(Key, Option<Value>)` instead.

  * vcf/record/genotypes/genotype/field: Remove `Field`.

    Use `(Key, Option<Value>)` instead.


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
