# Changelog

## Unreleased

### Added

  * vcf/record: Implement `TryFrom<&[u8]>`.

    This is a convenience implementation to allow constructing a single
    `vcf::Record` from a byte string.

## 0.61.0 - 2024-07-14

### Added

  * vcf/header/record/value/map/format/number: Add VCF 4.5 format numbers.

    VCF 4.5 introduces counts for local alternate bases (`LA`), local
    reference-alternate bases (`LR`), local samples (`LG`), ploidy (`P`), and
    base modifications (`M`).

  * vcf/header/record/value/map/format/definition: Add VCF 4.5 format
    definitions.

  * vcf/variant/record/samples/keys/key: Add VCF 4.5 format keys.

### Changed

  * vcf/header/record/value/map/format/definition/v4_4: Change phase set
    definitions to use number `P`.

    This was updated in [samtools/hts-specs@960c5fe].

[samtools/hts-specs@960c5fe]: https://github.com/samtools/hts-specs@960c5fe8cc6e3715f1321fb4de82f4d3824d77f6

## 0.60.0 - 2024-06-17

### Changed

  * vcf: Sync dependencies.

## 0.59.0 - 2024-06-06

### Changed

  * vcf/header/parser/record/value: Parse other header record values that start
    with a map prefix (`<`) as a string if there is no map identifier (i.e.,
    `ID=`) when the input is VCF < 4.3 ([#241]).

    This behavior is undefined in VCF < 4.3 and explicitly invalid in VCF >=
    4.3.

[#241]: https://github.com/zaeleus/noodles/issues/241

## 0.58.0 - 2024-05-31

### Changed

  * vcf/async/io/reader/header: Parse header line by line.

    The async header parser can now build a `vcf::Header` by parsing a raw
    header line by line. This makes it so that it is no longer required to read
    the entire raw header into memory before parsing.

  * vcf/io/writer/header/record/value: Disallow unstructured header records to
    be empty or start with a structured header record marker (`<`).

    See _The Variant Call Format Specification: VCFv4.4 and BCFv2.2_
    (2024-04-20) § 1.4 "Meta-information lines".

  * vcf/io/writer/record/reference_bases: Resolve IUPAC ambiguity codes
    ([#268]).

    Any occurrence of an IUPAC ambiguity code is resolved to its first single
    base specification. See _The Variant Call Format Specification: VCFv4.4 and
    BCFv2.2_ (2024-04-20) § 1.6.1.4 "Fixed fields: REF".

[#268]: https://github.com/zaeleus/noodles/issues/268

## 0.57.0 - 2024-05-16

### Added

  * vcf/header/record/value/map/format: Add `Number`.

    This is now separated from the `Number` that info records use.

### Changed

  * vcf/header: Move `Number` to `crate::header::record::value::map::info`
    module.

    This is now only used for info records.

## 0.56.0 - 2024-05-08

### Added

  * vcf/record/samples: Add lookup by sample name (`Samples::get`) and index
    (`Samples::get_index`).

  * vcf/variant/record_buf/samples: Add conversion to inner fields (`(Keys,
    Vec<Vec<Option<Value>>>)`).

  * vcf/variant/record_buf/samples: Add lookup by sample name (`Samples::get`).

### Changed

  * vcf/header/number: Expand variant names.

    Variants are renamed to descriptive names, i.e., `Number::A` =>
    `Number::AlternateBases`, `Number::R` => `Number::ReferenceAlternateBases`,
    and `Number::G` => `Number::Samples`. The single letter values are aliased
    to the long names.

### Fixed

  * vcf/io/writer/header/record: Write newlines for records in other
    collections ([#259]).

    This only affected collections with more than one record.

[#259]: https://github.com/zaeleus/noodles/issues/259

## 0.55.0 - 2024-05-02

### Removed

  * vcf/io/reader: Remove `Reader::virtual_position` and `Reader::seek`.

    Use the inner reader instead.

## 0.54.0 - 2024-04-22

### Added

  * vcf/variant/record/samples: Add series selector (`Samples::select`).

### Removed

  * vcf/header/file_format: Remove parser.

    An internal parser is used when parsing fileformat header records.

## 0.53.0 - 2024-04-11

### Changed

  * vcf/io/writer/record: Add validation for IDs, reference bases, alternate
    bases, filters, and info and samples keys.

  * vcf/io/writer/record: Percent-encode special characters in info and sample
    string values.

## 0.52.0 - 2024-04-04

### Added

  * vcf/header: Add string maps (`StringMaps`).

    This is moved from noodles-bcf.

  * vcf/io/indexed_reader: Add `IndexedReader::records`.

  * vcf/variant: Add `Record` trait to represent an opaque variant record.

    The variant record buffer is renamed to `RecordBuf`. This also introduces
    traits for field buffers.

  * vcf/variant/record_buf/samples/series/value: Add `Genotype` value.

    When decoding, samples genotype values (`GT`) can now be specialized as a
    list of allele position-phasing pairs.

### Changed

  * vcf: Move `VariantReader` and `VariantWriter` to `variant::io::Read` and
    `variant::io::Write`, respectively.

  * vcf: Move readers (`Reader` and `IndexedReader`) and writer (`Writer`) to
    `io` module.

  * vcf: Move `Record` to `variant::RecordBuf`.

  * vcf: Move `lazy::Record` to `Record`.

  * vcf/io/reader: Rename record to record buf and lazy record to record.

    This changes the following:

      * `Reader::read_record` => `Reader::read_record_buf`,
      * `Reader::records` => `Reader::record_bufs`,
      * `Records` => `RecordBufs`, and
      * `Reader::read_lazy_record` => `Reader::read_record`.

  * vcf/io/reader: Return an iterator over `vcf::Record` for queries
    (`Reader::query`).

  * vcf/variant/record: Rename position and end to variant start and variant
    end, respectively.

  * vcf/variant/io/write: Change `Write::write_record` to accept `&dyn
    crate::variant::Record`.

  * vcf/variant/record/info/field/key: Replace `Key` with a string.

  * vcf/variant/record/samples/keys/key: Replace `Key` with a string.

  * vcf/variant/record_buf: Rename chromosome to reference sequence name and
    genotypes to samples.

  * vcf/variant/record_buf: Simplify fields.

    This changes the following:

      * IDs => set of strings,
      * reference bases => string,
      * alternate bases => list of strings,
      * filters => set of strings, and
      * samples keys => set of strings.

  * vcf/variant/record_buf/builder: Remove validation on build
    (`Builder::build`).

    The record builder no longer checks whether a valid reference sequence
    name, position, or references bases is set.

  * vcf/variant/record_buf/samples: Samples can now be accessed by row or
    column.

    Samples are now represented like a data frame. They can be accessed per
    sample or per series.

### Removed

  * vcf/variant/record_buf: Remove `Position`.

    Record positions are replaced with `Option<core::Position>`, where
    `None` represents telomere start.

  * vcf/variant/record_buf: Remove `QualityScores`.

    Quality scores are replaced with `f32`.

  * vcf/variant/record_buf: Remove `RecordBuf::end`.

    Use `Record::end` instead.

  * vcf/variant/record_buf/samples: Remove `Samples::genotypes`.

    Use `Samples::select(key::GENOTYPE)` to get the raw values instead.

  * vcf/variant/record_buf/samples/sample/value/genotype/allele: Remove
    `Phasing`.

    Use `crate::variant::record::samples::series::value::genotype::Phasing`
    instead.

## 0.51.0 - 2024-03-28

### Changed

  * vcf: Sync dependencies.

## 0.50.0 - 2024-03-12

### Changed

  * vcf: Sync dependencies.

## 0.49.0 - 2024-01-25

### Changed

  * vcf: Sync dependencies.

### Fixed

  * vcf/lazy/record/bounds: Fix range for info field ([#223]).

  * vcf/reader/lazy_record: Disallow newlines to appear in fields ([#224]).

[#223]: https://github.com/zaeleus/noodles/pull/223
[#224]: https://github.com/zaeleus/noodles/pull/224

## 0.48.0 - 2023-12-14

### Changed

  * vcf: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * vcf/async/reader: Accept `csi::BinningIndex` for querying.

  * vcf/reader: Accept `csi::BinningIndex` for querying.

## 0.47.0 - 2023-11-14

### Changed

  * vcf: Sync dependencies.

## 0.46.0 - 2023-11-13

### Added

  * vcf: Add indexer (`vcf::index`) ([#214]).

    This is a convenience function to index a bgzip-compressed VCF file.

  * vcf/header/record/value/map: Add mutable getter for other fields
    (`Map::<I>::other_fields_mut`).

[#214]: https://github.com/zaeleus/noodles/issues/214

## 0.45.0 - 2023-11-02

### Changed

  * vcf/header/parser: Return error on duplicate structured line IDs.

## 0.44.0 - 2023-10-26

### Added

  * vcf/async/writer: Add a `shutdown` delegate.

    When the inner writer is buffered, a call to `AsyncWriter::shutdown` is
    required prior to drop.

  * vcf/header/parser: Add partial parser.

    This allows headers to be parsed line by line.

  * vcf/reader/builder: Add a compression method setter
    (`Builder::set_compression_method`).

    This allows the compression method to be overridden using
    `vcf::io::CompressionMethod`, instead of solely relying on the path
    extension.

    Change instantiations of `vcf::reader::Builder` to
    `vcf::reader::Builder::default()`.

  * vcf/writer/builder: Add `Builder::build_from_writer`.

  * vcf/writer/builder: Add a compression method setter.

    This allows the compression method to be overridden using
    `vcf::io::CompressionMethod`, instead of solely relying on the path
    extension.

    Change instantiations of `vcf::writer::Builder` to
    `vcf::writer::Builder::default()`.

### Changed

  * vcf/reader/builder: Change `Builder::build_from_reader` to accept
    implementations of `Read`.

  * vcf/reader/header: Parse header line by line.

    The header parser can now build a `vcf::Header` by parsing a raw header
    line by line. This makes it so that it is no longer required to read the
    entire raw header into memory before parsing.

## 0.43.0 - 2023-10-23

### Added

  * vcf/async: Add common methods to access the underlying I/O:
    `AsyncReader::get_ref`, `AsyncReader::get_mut`, `AsyncReader::into_inner`,
    `AsyncWriter::get_ref`, `AsyncWriter::get_mut`, and
    `AsyncWriter::into_inner`.

## 0.42.0 - 2023-10-19

### Changed

  * vcf: Sync dependencies.

## 0.41.0 - 2023-10-12

### Changed

  * vcf: Sync dependencies.

## 0.40.0 - 2023-09-21

### Changed

  * vcf/header/fmt: Add specialized META record serializer ([#203]).

    This no longer quotes the `Number`, `Type`, and `Values` field values.

  * vcf/header/parser/record/value/map/other: Add specialized PEDIGREE value
    parser ([#201]).

    When the input is VCF 4.2, this allows the `Child` or `Derived` field to
    act as the record ID in the value collection.

  * vcf/header/parser/record/value/map/other: Only parse `Values` field in VCF
    4.3+ META records.

[#201]: https://github.com/zaeleus/noodles/issues/201
[#203]: https://github.com/zaeleus/noodles/issues/203

## 0.39.0 - 2023-09-14

### Added

  * vcf/lazy/record: Implement `Clone` + `Debug` + `Eq` + `PartialEq`.

  * vcf/lazy/record: Add filters (`lazy::record::Filters`) and IDs
    (`lazy::record::Ids`) wrappers.

### Changed

  * vcf/lazy/record: Check missing field for quality score.

## 0.38.0 - 2023-08-31

### Added

  * vcf/async/reader: Add lazy record reader (`AsyncReader::read_lazy_record`;
    [#199]).

  * vcf/writer: Add builder (`vcf::writer::Builder`).

[#199]: https://github.com/zaeleus/noodles/pull/199

### Fixed

  * vcf/reader: Clear lazy record buffer before read.

  * vcf/record/genotypes/keys/key: Hash inner key.

  * vcf/record/info/field/key: Hash inner key.

## 0.37.0 - 2023-08-24

### Added

  * vcf/header/record/value/map/tag/other: Implement `PartialEq<&str>` for
    `Other<S>`.

  * vcf/indexed_reader: Add lazy record reader
    (`IndexedReader::read_lazy_record`).

  * vcf/lazy/record: Add info wrapper (`lazy::record::Info`).

### Changed

  * vcf/record/genotypes/keys/key: Allow full alphabet when the file format is
    < VCF 4.3.

  * vcf/record/info/field/key: Allow full alphabet when the file format is <
    VCF 4.3.

## 0.36.0 - 2023-08-17

### Added

  * vcf/header/record/value/map/tag: Implement `From<&str>` for `Tag<S>`.

  * vcf/header/record/value/map/contig/builder: Add URL setter
    (`Builder::set_url`).

  * vcf/lazy: Add a lazy record (`lazy::Record`).

    Lazy records are variant records that are lazily-evaluated. Their fields
    are not necessarily valid, but the buffer is guaranteed to be record-like.

  * vcf/reader: Add `Reader::read_lazy_record` to read lazy records.

  * vcf/record/position: Implement `PartialEq<core::Position>` and
    `PartialOrd<core::Position>` for `Position` ([#191]).

  * vcf/record/position: Implement `From<core::Position>` for `Position`.

[#191]: https://github.com/zaeleus/noodles/issues/191

### Changed

  * vcf/header/parser: Relax field order and value quoting rules.

    Required fields no longer have to be strictly ordered. Values can either be
    a raw or quoted string. Quoted strings are allowed to have some escape
    sequences.

  * vcf/header/record/value/map/contig/builder: Accept `Into<String>` for MD5
    checksum values (`Builder::set_md5`).

### Removed

  * vcf: Remove nom parser.

  * vcf/header: Remove `assembly`, `META`, and `pedigreeDB`.

    Use other records instead.

  * vcf/header/record/value/map: Remove `TryFrom<Vec<String, String>>`.

    Use the map builders instead.

## 0.35.0 - 2023-08-03

### Changed

  * vcf: Sync dependencies.

## 0.34.0 - 2023-07-06

### Changed

  * vcf: Update to indexmap 2.0.0.

  * vcf/record/info: Change `Info::get_index_mut` to return a reference to the
    key.

    This was previously a _mutable_ reference to the key, but this is only safe
    to change if the new value matches the previous key's hash and quality,
    which is (very likely) never the case for info keys.

## 0.33.0 - 2023-06-29

### Changed

  * vcf/header/file_format: Make methods `const`.

## 0.32.0 - 2023-06-15

### Changed

  * vcf: Sync dependencies.

## 0.31.0 - 2023-06-01

### Changed

  * vcf: Sync dependencies.

## 0.30.0 - 2023-05-18

### Changed

  * vcf: Sync dependencies.

## 0.29.0 - 2023-05-04

### Added

  * vcf/record/genotypes/sample/value: Add conversions from raw types.

  * vcf/record/info/field/value: Add conversions from raw types.

### Changed

  * vcf/header/file_format: The default file format is now VCF 4.4.

  * vcf/record/genotypes/sample/value: Add array value wrapper (`Array`).

    This replaces `Value::*Array(_)` with `Value::Array(Array)`.

  * vcf/record/info/field/value: Add array value wrapper (`Array`).

    This replaces `Value::*Array(_)` with `Value::Array(Array)`.

## 0.28.0 - 2023-04-27

### Added

  * vcf/header/record/key/other: Implement `FromStr` for `Other`.

  * vcf/header/record/value/map/contig: Add getter (`Map::<Contig>::url`) and
    mutable getter (`Map::<Contig>::url_mut`) for URL.

  * vcf/header/record/value/map/tag/other: Implement `FromStr` for `Other`.

  * vcf/record/genotypes/sample: Implement `Debug` and `PartialEq`
    for `Sample`.

  * vcf/record/info: Change getters to accept keys that implement
    `Equivalent<Key>`.

### Changed

  * vcf/header: Change other records to a collection of either unstructured or
    structured records.

    This changes `OtherRecords` to a be `IndexMap<key::Other,
    value::Collection>`, where the collection can be a list of strings or a map
    of `Map<Other>` values. Change usages of `header::record::value::Other` to
    `header::record::Value`.

  * vcf/header: `Header::insert` is now fallible.

    This can fail when the value type does not match the collection type for
    the given key.

  * vcf/header/format: Move `Key` under `record::genotypes::keys` module.

  * vcf/header/info: Move `Key` under `record::info::field` module.

  * vcf/header/number: Change default to `Number::Count(1)`.

    The default type definition is simplified to a single string.

  * vcf/header/record/value/map: Increase the visibility of `Tag` for inner
    maps.

  * vcf/header/record/value/map/builder: Accept `Into<String>` for other field
    values (`Builder::insert`).

  * vcf/header/record/value/map/format: Implement `From<(FileFormat, &Key)>`.

  * vcf/header/record/value/map/format/ty: Implement `Default`.

    `Type::String` is now the default format type.

  * vcf/header/record/value/map/info: Implement `From<(FileFormat, &Key)>`.

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

  * vcf/header/record/value/map: Remove `TryFromFieldsError`.

    Use the inner record `ParseError` instead, e.g., `map::info::ParseError`,
    etc.

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
