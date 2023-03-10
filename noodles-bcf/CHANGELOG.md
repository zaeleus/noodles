# Changelog

## Unreleased

### Added

  * bcf/writer: Implement `VariantWriter` ([#150]).

[#150]: https://github.com/zaeleus/noodles/issues/150

### Fixed

  * bcf/writer/vcf_record/genotypes: Fix writing dropped values ([#151]).

    Trailing field values are allowed to be missing.

[#151]: https://github.com/zaeleus/noodles/issues/151

## 0.21.0 - 2023-03-03

### Added

  * bcf/reader: Implement `vcf::VariantReader` ([#149]).

[#149]: https://github.com/zaeleus/noodles/pull/149

## 0.20.0 - 2023-02-03

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.64.0.

  * bcf/record/info: Add `Info::iter` for an iterator over
    `Key`-`Option<Value>` pairs.

  * bcf/record/info: Change `Info::values` to return an iterator over
    `Option<Value>` values.

    Use `Info::iter` for fields.

  * bcf/record/info: Change `Info::get` to return an `Option<Value>` value
    instead of a field.

### Removed

  * bcf/header: Remove `StringMap`.

    This was deprecated in noodles-bcf 0.11.0. Use
    `noodles_bcf::header::string_maps::StringMap` instead.

## 0.19.2 - 2022-12-02

### Fixed

  * bcf/writer/vcf_writer/genotypes: Add FORMAT GT encoder for missing alleles
    (i.e., `.`).

  * bcf/writer/vcf_writer/genotypes: Pad integer vector for genotypes.

## 0.19.1 - 2022-11-29

### Fixed

  * bcf/header/string_maps: Use file format as context when parsing records
    (#137).

  * bcf/writer/vcf_writer/genotypes: Add FORMAT GT encoder ([#135]).

[#135]: https://github.com/zaeleus/noodles/issues/135
[#137]: https://github.com/zaeleus/noodles/issues/137

## 0.19.0 - 2022-11-18

### Changed

  * bcf: Sync dependencies.

## 0.18.0 - 2022-10-28

### Changed

  * bcf: Sync dependencies.

## 0.17.0 - 2022-10-20

### Changed

  * bcf: Sync dependencies.

## 0.16.0 - 2022-09-29

### Changed

  * bcf: Sync dependencies.

## 0.15.0 - 2022-08-16

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.57.0.

### Removed

  * bcf/async/reader: Remove `Builder` and `Reader::builder`.

    Use `bgzf::async::reader::Builder` and construct with
    `bcf::async::Reader::from` instead.

## 0.14.0 - 2022-07-05

### Changed

  * bcf: Sync dependencies.

## 0.13.3 - 2022-06-08

### Fixed

  * bcf: Sync dependencies.

## 0.13.2 - 2022-03-29

### Fixed

  * bcf/async/reader/record: Read site from buffer ([#79]).

[#79]: https://github.com/zaeleus/noodles/issues/79

## 0.13.1 - 2022-03-02

### Fixed

  * bcf: Sync dependencies.

## 0.13.0 - 2022-02-17

### Added

  * bcf: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

  * bcf/async/reader: Add query stream.

  * bcf/record: Add type alias for chromosome ID: `ChromosomeId`.

### Fixed

  * bcf/writer/vcf_record/site: Avoid overflow when the alternate allele count
    is `usize::MAX`.

## 0.12.0 - 2022-01-27

### Added

  * bcf/async/reader: Add common methods to access the underlying reader:
    `get_ref`, `get_mut`, and `into_inner`.

  * bcf/reader/record/info: Handle reading INFO array values with missing raw
    values.

  * bcf/writer/vcf_record/site/info: Handle writing INFO array values with
    missing raw values.

    This allows reading/writing INFO array values that have missing raw values,
    e.g., `AC=8,.`.

### Changed

  * bcf/record: Change chromosome ID to a `usize`.

    `CHROM` is used as an index.

  * bcf/reader: Use contig string map to find contig index when querying.

    The chromosome ID maps to a contig string map entry, not the header
    contigs.

  * bcf/writer/vcf_record/site: Ensure the record position is <= its end
    position when writing `rlen`.

### Fixed

  * bcf/reader/record/info: Split on delimiter when reading a character array.

  * bcf/writer/record: Write filters when present.

  * bcf/writer/vcf_record/site: Skip writing alternate bases when empty.

  * bcf/writer/vcf_record/site: Filters with indices larger than 127 are now
    valid.

  * bcf/writer/vcf_record/site: Fix `rlen` calculation.

## 0.11.0 - 2022-01-13

### Added

  * bcf/async/reader: Add conversion from `R` into `Reader<R>`.

### Changed

  * bcf/header: Split `StringMap` from `StringMaps` (#64).

    Use `bcf::header::string_maps::StringStringMap` as a replacement to what
    `StringMap` was used as. It can typically be accessed using
    `StringMaps::strings`.

  * bcf/header/string_maps: Parsing can now fail with
    `vcf::header::ParseError::StringMapPositionMismatch` if the string map
    position of an entry and record-defined IDX field value do not match.

  * bcf/header/string_maps: If present, the IDX field is used to determine the
    position of the entry in the string map ([#64]).

  * bcf/header/string_maps: Parsing (`FromStr`) and conversion
    (`From<&vcf::Header>`) now builds a contig string map (also known as a
    "dictionary of contigs") (#64).

    It can be accessed using `StringMaps::contigs`.

  * bcf/header/string_maps/string_map: The default implementation now creates
    an empty map.

    This used to create a default map with a `PASS` entry at position 0.

  * bcf/record: `Record::try_into_vcf_record` now takes `StringMaps` instead of
    `StringMap` (#64).

  * bcf/record: The following now take `StringStringMap` instead of
    `StringMap` (#64):

    * `Filters::try_into_vcf_record_filters`
    * `Info::try_into_vcf_record_info`
    * `Info::get`
    * `Info::values`
    * `Info::try_into_vcf_record_genotypes`

  * bcf/record/convert: Use contig string map to find chromosome name (#64).

  * bcf/writer: `Writer::write_vcf_record` now takes `StringMaps` instead of
    `StringMap`.

  * bcf/writer/vcf_record: Use contig string map to find index of contig (#64).

    This previously used the position of the contig in the header, but BCF
    allows records to override their positions using the `IDX` field, requiring
    a contig string map.

[#64]: https://github.com/zaeleus/noodles/issues/64

### Removed

  * bcf/header/string_maps/string_map: Remove `Deref<Target =
    IndexSet<String>>` for `StringMap`.

    `StringMap` is no longer backed by an `IndexMap`.

## 0.10.0 - 2021-12-16

### Added

  * bcf/reader: Add common methods to access the underlying reader: `get_ref`,
    `get_mut` and `into_inner`.

  * bcf/writer: Add record writer (`Writer::write_record`).

### Fixed

  * bcf/writer/vcf_record/site: Fix allele count calculation.

    This overcounted by 1 when there were no alternate bases.

## 0.9.0 - 2021-12-02

### Added

  * bcf/record/genotypes: Add method to clear all fields (`Genotypes::clear`).

  * bcf/record/info: Add methods to wrap an INFO buffer (`Info::new`), clear
    all fields (`Info::clear`), retrieve a field by key (`Info::get`; [#52]),
    and iterate all fields (`Info::values`).

  * bcf/reader: Add conversion from `R` into `Reader<R>`.

  * bcf/writer: Add conversion from `W` into `Writer<W>`.

  * bcf/writer: Add common methods to access the underlying writer: `get_mut`
    and `into_inner`.

[#52]: https://github.com/zaeleus/noodles/issues/52

### Changed

  * bcf/header/string_map: Disable type checking when the file format is < VCF
    4.3.

  * bcf/reader/record/genotypes: Always use key from header.

  * bcf/reader/record/site/info: Always use key from header.

  * bcf/record/genotypes: Return `Genotypes` from
    `Genotypes::try_into_vcf_record_genotypes`.

    Use `vcf::record::Genotypes::keys` for the genotypes keys.

### Fixed

  * bcf/reader/record/site/info: Allow value to be optional (#56).

[#56]: https://github.com/zaeleus/noodles/issues/56

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
