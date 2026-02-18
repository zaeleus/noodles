# Changelog

## 0.82.0 - 2026-02-18

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.88.0.

## 0.81.0 - 2025-12-18

### Changed

  * bcf/io/reader/query: Change `Query` to a reader.

    `Query` is now a reader adapter, rather than an iterator. Use
    `Query::records` for an iterator.

## 0.80.0 - 2025-12-11

### Removed

  * bcf: Remove deprecated items.

    The following items are removed:

      * `AsyncReader` (deprecated since 0.76.0; use `bcf::r#async::io::Reader`
        instead) and
      * `AsyncWriter` (0.76.0; `bcf::r#async::io::Writer`).

## 0.79.0 - 2025-11-13

### Changed

  * bcf: Sync dependencies.

## 0.78.0 - 2025-08-25

### Changed

  * bcf: Sync dependencies.

## 0.77.0 - 2025-07-12

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.85.0.

  * bcf/record/encoder/site: Use variant span for `rlen` field.

    Writing would previously fail for telomeric breakend records (i.e., records
    with missing positions).

## 0.76.0 - 2025-05-29

### Deprecated

  * bcf: Deprecate async re-exports (`AsyncReader` and
    `AsyncWriter`).

    Use `bcf::r#async::io::Reader` and `bcf::r#async::io::Writer`
    instead.

## 0.75.0 - 2025-05-16

### Changed

  * bcf: Sync dependencies.

## 0.74.0 - 2025-04-13

### Changed

  * bcf: Sync dependencies.

## 0.73.0 - 2025-04-06

### Changed

  * bcf: Sync dependencies.

## 0.72.0 - 2025-03-08

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.81.0.

### Fixed

  * bcf/record/filters: Fix indices iterator.

    This fixes decoding lists with indices > `i8::MAX`.

## 0.71.0 - 2025-02-17

### Added

  * bcf/fs: Add indexer (`index`).

### Fixed

  * bcf/record/samples/series: Fix series length ([#328]).

    This incorrectly used the value size instead of the sample count.

[#328]: https://github.com/zaeleus/noodles/issues/328

## 0.70.0 - 2025-02-06

### Changed

  * bcf: Sync dependencies.

## 0.69.0 - 2025-01-24

### Changed

  * bcf: Sync dependencies.

## 0.68.0 - 2025-01-23

### Changed

  * bcf: Sync dependencies.

## 0.67.0 - 2025-01-19

### Added

  * bcf/io/reader: Add a header reader adapter (`header::Reader`).

    This can be used to read the raw fields of a BCF header, e.g., the raw VCF
    header.

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.73.0.

  * bcf/async/io/reader/header: Parse header line by line.

    The async header reader now builds a `vcf::Header` by parsing a raw header
    line by line. This makes it so that it is no longer required to read the
    entire raw header into memory before parsing.

  * bcf/async/io/reader/header: Discard VCF header padding.

## 0.66.0 - 2024-12-20

### Changed

  * bcf: Sync dependencies.

## 0.65.0 - 2024-12-12

### Changed

  * bcf: Sync dependencies.

## 0.64.0 - 2024-11-07

### Changed

  * bcf: Sync dependencies.

## 0.63.0 - 2024-10-22

### Changed

  * bcf/record/samples/series/value/genotype: Only infer first allele phasing
    for VCF < 4.4.

    See § 6.3.3.9 "Type encoding: Genotype (GT) field" (2024-06-28): "When
    processing VCF version 4.3 or earlier files, the phasing of the first
    allele should be treated as missing and inferred from the remaining
    values."

## 0.62.0 - 2024-09-26

### Changed

  * bcf: Sync dependencies.

## 0.61.0 - 2024-09-12

### Changed

  * bcf: Sync dependencies.

## 0.60.0 - 2024-09-09

### Changed

  * bcf: Sync dependencies.

## 0.59.1 - 2024-09-04

### Fixed

  * bcf/record/codec/decoder/position: Fix reading position at max position.

    This would previously overflow and error instead of returning a properly
    normalized `Position`.

## 0.59.0 - 2024-08-04

### Changed

  * bcf: Sync dependencies.

## 0.58.0 - 2024-07-14

### Changed

  * bcf: Sync dependencies.

## 0.57.0 - 2024-06-17

### Changed

  * bcf: Sync dependencies.

## 0.56.0 - 2024-06-06

### Changed

  * bcf: Sync dependencies.

## 0.55.0 - 2024-05-31

### Changed

  * bcf/record: Resolve info field and samples series values using header type
    definitions ([#266]).

[#266]: https://github.com/zaeleus/noodles/issues/266

## 0.54.0 - 2024-05-16

### Changed

  * bcf/record/samples: Add lookup by sample name (`Samples::get`).

## 0.53.0 - 2024-05-08

### Changed

  * bcf: Sync dependencies.

## 0.52.0 - 2024-05-02

### Removed

  * bcf/io/reader: Remove `Reader::virtual_position` and `Reader::seek`.

    Use the inner reader instead.

## 0.51.0 - 2024-04-22

### Added

  * bcf/record/samples: Implement series selector (`Samples::select`).

### Fixed

  * bcf/io/reader/record: Return EOF when site length signals EOF ([#255]).

[#255]: https://github.com/zaeleus/noodles/issues/255

## 0.50.0 - 2024-04-11

### Changed

  * bcf/io/reader/record: Disallow partial site lengths.

    A partial site length previously would return EOF. This is now more
    strict, requiring a complete site length or EOF.

## 0.49.0 - 2024-04-04

### Added

  * bcf/async/io: Add async writer (`r#async::io::Writer`).

  * bcf/record: Add reference sequence name getter
    (`Record::reference_sequence_name`).

  * bcf/record: Add `Samples` wrapper.

    This allows reading of the samples without converting to
    `vcf::variant::record_buf::Samples`.

### Changed

  * bcf: Move readers (`Reader` and `IndexedReader`) to writer (`Writer`) to
    `io` module.

  * bcf: Move lazy record to record.

    Use `bcf::Record` instead of `bcf::lazy::Record`.

  * bcf/async/io/reader: Change `Reader::read_header` to return a parsed header
    (`vcf::Header`).

    This no longer returns a raw string. String maps are also parsed and
    attached to the reader. Use `Reader::string_maps` to get a reference to
    the string maps.

  * bcf/io/reader: Rename record to record buf.

    This changes the following:

      * `Reader::read_record` => `Reader::read_record_buf`,
      * `Reader::records` => `Reader::record_bufs`,
      * `Reader::read_lazy_record` => `Reader::read_record`, and
      * `Reader::lazy_records` => `Reader::records`.

  * bcf/io/reader: Return an iterator over `bcf::Record` for queries
    (`Reader::query`).

  * bcf/record: Rename chromosome to reference sequence name, position to
    variant start, and genotypes to samples.

  * bcf/record :Increase the visibilities of the bases methods
    (`Record::reference_bases` and `Record::alternate_bases`).

### Removed

  * bcf/async/io/reader: Remove `Reader::read_file_format`.

    `Reader::read_header` now includes checking the magic number and reading
    the file format version.

  * bcf/header: Remove `StringMaps`.

    This is moved to noodles-vcf.

  * bcf/record: Remove `Record::try_into_vcf_record`.

    Use `vcf::variant::RecordBuf::try_from_variant_record` instead.

## 0.48.0 - 2024-03-28

### Changed

  * bcf: Sync dependencies.

## 0.47.0 - 2024-03-12

### Changed

  * bcf: Sync dependencies.

## 0.46.0 - 2024-01-25

### Changed

  * bcf: Sync dependencies.

## 0.45.0 - 2023-12-14

### Added

  * bcf/indexed_reader: Add index getter (`IndexedReader::index`).

  * bcf/reader/builder: Add `Builder::build_from_reader`.

  * bcf/reader/builder: Add a compression method setter
    (`Builder::set_compression_method`).

    This allows the compression method to be overridden using
    `bcf::io::CompressionMethod`, instead of assuming the input is
    always bgzip-compressed.

    Change instantiations of `bcf::reader::Builder` to
    `bcf::reader::Builder::default()`.

  * bcf/writer: Add builder (`bcf::writer::Builder`).

### Changed

  * bcf: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * bcf: Increase the visibility of the `writer` module.

  * bcf/async/reader: Accept `csi::BinningIndex` for querying.

  * bcf/reader: Accept `csi::BinningIndex` for querying.

## 0.44.0 - 2023-11-14

### Changed

  * bcf: Sync dependencies.

## 0.43.0 - 2023-11-13

### Changed

  * bcf: Sync dependencies.

## 0.42.0 - 2023-11-02

### Changed

  * bcf/reader/header: Parse header line by line.

    The header parser can now build a `vcf::Header` and `StringMaps` by parsing
    a raw header line by line. This makes it so that it is no longer required
    to read the entire raw header into memory before parsing.

## 0.41.0 - 2023-10-26

### Changed

  * bcf: Sync dependencies.

## 0.40.0 - 2023-10-23

### Added

  * bcf/reader: Add builder (`bcf::reader::Builder`).

## 0.39.0 - 2023-10-19

### Changed

  * bcf: Sync dependencies.

## 0.38.0 - 2023-10-12

### Changed

  * bcf: Sync dependencies.

## 0.37.0 - 2023-09-21

### Changed

  * bcf: Sync dependencies.

## 0.36.0 - 2023-09-14

### Removed

  * bcf/lazy/record: Remove mutable reference getters.

    Lazy records were not meant to be mutable. This removes
    `lazy::Record::position_mut`, `lazy::Record::quality_score_mut`, and
    `lazy::Record::ids_mut`.

## 0.35.0 - 2023-08-31

### Changed

  * bcf: Sync dependencies.

## 0.34.0 - 2023-08-24

### Changed

  * bcf: Sync dependencies.

## 0.33.0 - 2023-08-17

### Changed

  * bcf: Sync dependencies.

## 0.32.0 - 2023-08-03

### Added

  * bcf/indexed_reader/builder: Add index setter (`Builder::set_index`).

  * bcf/indexed_reader/builder: Add `Builder::build_from_reader`.

## 0.31.0 - 2023-07-06

### Changed

  * bcf: Update to indexmap 2.0.0.

### Fixed

  * bcf/reader/query: Fix using parsed string maps for query ([#181]).

[#181]: https://github.com/zaeleus/noodles/issues/181

## 0.30.0 - 2023-06-29

### Added

  * bcf/indexed_reader: Add builder (`indexed_reader::Builder`).

## 0.29.0 - 2023-06-15

### Added

  * bcf: Add an indexed reader (`IndexedReader`).

## 0.28.0 - 2023-06-01

### Changed

  * bcf: Sync dependencies.

## 0.27.0 - 2023-05-18

### Changed

  * bcf: Sync dependencies.

## 0.26.0 - 2023-05-11

### Removed

  * bcf/reader: Remove `Reader::read_file_format`.

    `Reader::read_header` now includes checking the magic number and reading
    the file format version. Calling `Reader::read_file_format` is no longer
    necessary.

  * bcf/writer: Remove `Writer::write_file_format`.

    `Writer::write_header` now includes writing the magic number and file
    format version. Calling `Writer::write_file_format` is no longer necessary.

## 0.25.0 - 2023-05-04

### Changed

  * bcf: Sync dependencies.

## 0.24.0 - 2023-04-27

### Changed

  * bcf: Rename `Record` to `lazy::Record`.

  * bcf/async/reader: Rename `Reader::read_record` to
    `Reader::read_lazy_record`.

  * bcf/async/reader: Rename `Reader::records` to `Reader::lazy_records`.

  * bcf/header/string_maps: Use the `IDX` field, if present, when building
    string maps from a `vcf::Header`.

    This now makes conversion fallible. Use `TryFrom<&vcf::Header>` instead.

  * bcf/reader: Change `Reader::read_header` to return a parsed header
    (`Header`).

    This no longer returns a raw string. `StringMaps` are also parsed and
    attached to the reader. Use `Reader::string_maps` to get a reference to the
    string maps.

  * bcf/reader: Change `Reader:read_record` to read a parsed record.

    This now requires the header and strings map. If reading a lazy record, use
    `Reader::read_lazy_record` instead.

  * bcf/reader: Rename `Reader::records` to `Reader::lazy_records`.

  * bcf/reader/query: Return `vcf::Record` instead of `bcf::lazy::Record`.

    This now requires `Reader::query` to receive a reference to a
    `vcf::Header`.

  * bcf/writer: Build string maps upon writing header.

    `StringMaps` are now built from the header upon calling
    `Writer::write_header`. This is subsequently used as context when writing
    records.

  * bcf/writer: Rename `Writer::write_vcf_record` to `Writer::write_record`.

### Removed

  * bcf/writer: Remove `Writer::write_record(lazy::Record)`.

    Writing lazy records is not supported.

## 0.23.0 - 2023-04-06

### Changed

  * bcf: Sync dependencies.

## 0.22.0 - 2023-03-14

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
