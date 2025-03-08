# Changelog

## 0.50.0 - 2025-03-08

### Changed

  * fasta: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.49.0 - 2025-02-06

### Changed

  * fasta: Sync dependencies.

## 0.48.0 - 2025-01-24

### Changed

  * fasta/io: Move file indexer (`index`) to top-level `fs` module.

### Deprecated

  * fasta/io: Deprecate `index`.

    Use `fasta::fs::index` instead.

## 0.47.0 - 2025-01-19

### Added

  * fasta/fai/async/fs: Add convenience read function
    (`fai::r#async::fs::read`).

  * fasta/fai/async/io: Add writer (`Writer`).

  * fasta/fai/fs: Add convenience write function (`fai::fs::write`).

### Changed

  * fasta: Raise minimum supported Rust version (MSRV) to 1.73.0.

  * fasta/fai: Move convenience `read` function to `fs` module.

### Deprecated

  * fasta/fai: Deprecate `read`.

    Use `fai::fs::read` instead.

## 0.46.0 - 2024-12-12

### Changed

  * fasta: Sync dependencies.

## 0.45.0 - 2024-10-22

### Added

  * fasta/async/io: Add async writer (`fasta::r#async::io::Writer`).

  * fasta/io/writer/builder: Add build from path (`Builder::build_from_path`).

### Deprecated

  * fasta: Deprecate async re-export (`AsyncReader`).

    Use `fasta::r#async::io::Reader` instead.

## 0.44.0 - 2024-09-26

### Changed

  * fasta/fai: Move reader (`Reader`) and writer (`Writer`) to `io` module.

  * fasta/fai/async: Move reader (`Reader`) to `io` module.

### Deprecated

  * fasta/fai: Deprecate `Reader` and `Writer`.

    Use `fai::io::Reader` and `fai::io::Writer`, respectively, instead.

  * fasta/fai/async: Deprecate `Reader` and `Writer`.

    Use `fai::r#async::io::Reader` instead.

## 0.43.0 - 2024-09-04

### Changed

  * fasta/io/writer/builder: Rename `Builder::build_with_writer`
    to `Builder::build_from_writer`.

### Deprecated

  * fasta/io/writer/builder: Deprecate `Builder::build_with_writer`.

    Use `Builder::build_from_writer` instead.

## 0.42.0 - 2024-08-04

### Added

  * fasta/io/writer: Add additional methods to get the underlying writer
    (`Writer::get_mut` and `Writer::into_inner`).

## 0.41.0 - 2024-07-14

### Added

  * fasta/async/io/reader: Add common methods to access the underlying reader.

  * fasta/fai/record: Implement `Clone`.

  * fasta/sequence: Add `Record` trait.

### Changed

  * fasta/fai/index: Add wrapper (`Index`).

    Use `AsRef<Vec<Record>>` to get the list of records.

## 0.40.0 - 2024-06-17

### Added

  * fasta/fai: Add common methods to access the underlying I/O ([#269]).

[#269]: https://github.com/zaeleus/noodles/issues/269

## 0.39.0 - 2024-05-31

### Changed

  * fasta: Move readers (`Reader` and `IndexedReader`), writer (`Writer`), and
    indexer (`index`) to `io` module.

  * fasta/record/sequence/complement: Allow lowercase bases ([#267]).

[#267]: https://github.com/zaeleus/noodles/pull/267

### Deprecated

  * fasta: Deprecate `fasta::Reader`, `fasta::IndexedReader`, `fasta::Writer`,
    and `fasta::index`.

    These are moved under the `fasta::io` module.

## 0.38.0 - 2024-05-16

### Changed

  * fasta: Sync dependencies.

## 0.37.0 - 2024-05-08

### Changed

  * fasta: Sync dependencies.

## 0.36.0 - 2024-05-02

### Changed

  * fasta: Sync dependencies.

## 0.35.0 - 2024-03-28

### Changed

  * fasta: Sync dependencies.

## 0.34.0 - 2024-03-12

### Changed

  * fasta: Sync dependencies.

## 0.33.0 - 2024-02-22

### Changed

  * fasta/repository/adapter: Add `Send` and `Sync` constraints to `Adapter`.

### Fixed

  * fasta/repository: Make `Repository` thread-safe ([#234]).

  `Repository` was intended to allow concurrent access, which is why the
  adapter has a readers-writer lock.

[#234]: https://github.com/zaeleus/noodles/issues/234

## 0.32.0 - 2024-01-25

### Changed

  * fasta/record/definition: Change fields to byte strings.

## 0.31.0 - 2023-12-14

### Changed

  * fasta: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.30.0 - 2023-10-12

### Changed

  * fasta: Sync dependencies.

## 0.29.0 - 2023-08-31

### Changed

  * fasta: Sync dependencies.

## 0.28.0 - 2023-08-17

### Changed

  * fasta: Sync dependencies.

## 0.27.0 - 2023-07-27

### Changed

  * fasta/reader/sequence: Skip empty lines ([#189]).

[#189]: https://github.com/zaeleus/noodles/issues/189

## 0.26.0 - 2023-07-20

### Changed

  * fasta/reader: Increase visibility of `sequence::Reader`.

## 0.25.0 - 2023-06-15

### Changed

  * fasta: Sync dependencies.

## 0.24.0 - 2023-06-08

### Added

  * fasta/reader: Add a sequence reader (`fasta::reader::sequence::Reader`)
    ([#101]).

    This is created by calling `fasta::Reader::sequence_reader`. It is used for
    lower-level reading of the sequence.

[#101]: https://github.com/zaeleus/noodles/issues/101

## 0.23.0 - 2023-06-01

### Changed

  * fasta: Sync dependencies.

### Deprecated

  * fasta/fai/record: Deprecate `Record::len`.

    Use `Record::length` instead.

## 0.22.0 - 2023-05-04

### Changed

  * fasta/fai/record: Accept `Into<String>` for the record name.

## 0.21.0 - 2023-04-27

### Changed

  * fasta: Sync dependencies.

## 0.20.0 - 2023-03-14

### Added

  * fasta/indexed_reader: Add getter for index (`IndexedReader::index`).

## 0.19.0 - 2023-03-03

### Changed

  * fasta/reader: Improve performance of querying small regions ([#146]).

    The reader no longer needs to read the entire sequence when querying.

  * fasta/reader/builder: Add `bgz` as a known BGZF extension ([#144]).

[#144]: https://github.com/zaeleus/noodles/pull/144
[#146]: https://github.com/zaeleus/noodles/issues/146

## 0.18.0 - 2023-02-03

### Added

  * fasta/fai/record: Implement `std::error::Error::source` for `ParseError`.

### Changed

  * fasta: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * fasta/fai/record: Remove `Record::reference_sequence_name`.

    This was deprecated in noodles-fasta 0.3.0. Use `Record::name`
    instead.

  * fasta/record: Remove `Record::reference_sequence_name`.

    This was deprecated in noodles-fasta 0.2.0. Use `Record::name`
    instead.

  * fasta/record/definition: Remove
    `Definition::reference_sequence_name`.

    This was deprecated in noodles-fasta 0.3.0. Use `Definition::name`
    instead.

  * fasta/record/definition: Remove
    `ParseError::MissingReferenceSequenceName`.

    This was deprecated in noodles-fasta 0.6.0. Use
    `ParseError::MissingName` instead.


## 0.17.0 - 2022-11-18

### Changed

  * fasta: Sync dependencies.

## 0.16.0 - 2022-10-28

### Added

  * fasta/writer/builder: Implement `Default`.

### Changed

  * fasta/writer/builder: `Builder` no longer holds a writer.

### Removed

  * fasta/writer: Remove `Writer::builder`.

    Use `writer::Builder::default` instead.

  * fasta/writer/builder: Remove `Builder::build`.

    Use `Builder::build_with_writer` instead.

## 0.15.0 - 2022-10-20

### Changed

  * fasta: Sync dependencies.

### Unreleased

## Changed

  * fasta: Split indexed reader from reader.

    `reader::Builder` no longer attempts to load associated indices. This
    functionality is separated into `IndexedReader`, which now guarantees
    associated indices are loaded for querying.

    Changes usages of `fasta::reader::Builder` to
    `fasta::indexed_reader::Builder` if it is known querying is necessary.

## 0.14.0 - 2022-09-29

### Added

  * fasta/reader: Add common methods to access the underlying reader
    (`Reader::get_ref`, `Reader::get_mut`, and `Reader::into_inner`).

  * fasta/reader: Add builder (`fasta::reader::Builder`).

    The builder is able to construct a reader from a path
    (`Builder::build_from_path`), which can open raw FASTA files (`*.fa`) and
    bgzipped FASTA (`*.fa.gz`) files. If an associated index (`*.fai`) exists,
    it is loaded to allow querying.

### Changed

  * fasta/reader: `Reader::query` no longer takes a `fai::Index` as input.

    Use an indexed reader via `reader::Builder` to load or set an associated
    index instead.

### Removed

  * fasta/reader: Remove `seek` and `virtual_position` delegates.

    Use the inner reader instead.

  * fasta/repository/adapters/indexed_reader: Remove `Builder`.

    Use `fasta::reader::Builder` and construct `IndexedReader` with a
    `fasta::Reader` instead.

## 0.13.0 - 2022-08-16

### Changed

  * fasta: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.12.0 - 2022-07-05

### Changed

  * fasta: Sync dependencies.

## 0.11.0 - 2022-06-08

### Added

  * fasta/record/sequence: Implement `FromIterator<u8>`.

  * fasta/record/sequence: Add iterator to complement a sequence
    (`Sequence::complement`) ([#86]).

  * fasta/record/sequence: Add method to return a slice as a `Sequence`
    (`Sequence::slice`).

  * fasta/writer: Add builder ([#87]).

    This allows overriding the line base count.

[#86]: https://github.com/zaeleus/noodles/issues/86
[#87]: https://github.com/zaeleus/noodles/issues/87

### Changed

  * fasta/indexer: Sequence lines no longer strip end-of-line ASCII whitespace.

    This previously would take a line (i.e., characters up to `\n`) and strip
    [ASCII whitespace] from the end. To be consistent with the FASTA reader,
    all characters in the line are now considered bases. This leads to a ~20%
    performance improvement by avoiding having to copy the line buffer.

[ASCII whitespace]: https://infra.spec.whatwg.org/#ascii-whitespace

## 0.10.0 - 2022-04-14

### Added

  * fasta/repository/adapters/indexed_reader: Add convenience builder for
    `BufReader<File>`.

## 0.9.0 - 2022-03-29

### Added

  * fasta/record/sequence: Add indexing by `Position`.

  * fasta/record/sequence: Implement `From<Bytes>`.

### Changed

  * fasta/repository: Box the adapter.

    `fasta::Repository` no longer carries an adapter generic.

  * fasta/repository: Implement `Default`.

  * fasta/repository/adapters: Add an empty adapter.

    This may be useful to create a repository that is never used.

## 0.8.0 - 2022-03-02

### Added

  * fasta: Add async reader (`fasta::AsyncReader`).

  * fasta: Add sequence repository (`fasta::Repository`).

    A repository is a concurrent cache that uses a storage adapter to lookup
    and load sequence data.

  * fasta/fai: Add async reader (`fai::AsyncReader`).

### Changed

  * fasta/record/sequence: `Sequence` is now backed by a `Bytes` buffer.

    This allows for zero-copy cloning of the sequence or slices of the
    sequence.

## 0.7.0 - 2022-02-17

### Added

  * fasta: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.6.0 - 2022-01-27

### Deprecated

  * fasta/record/definition: Deprecate
    `ParseError::MissingReferenceSequenceName`.

    Use `ParseError::MissingName` instead.

## 0.5.2 - 2022-01-13

### Fixed

  * fasta: Sync dependencies.

## 0.5.1 - 2021-12-09

### Fixed

  * fasta: Sync dependencies.

## 0.5.0 - 2021-12-02

### Added

  * fasta/record/definition: Accept `Into<String>` for name.

### Changed

  * fasta/record: Wrap sequence.

    Use `Sequence::as_ref` to get the underlying list.

## 0.4.1 - 2021-11-18

### Fixed

  * fasta: Sync dependencies.

## 0.4.0 - 2021-11-11

### Changed

  * fasta: Update to Rust 2021.

## 0.3.1 - 2021-10-16

### Fixed

  * fasta: Sync dependencies.

## 0.3.0 - 2021-10-01

### Deprecated

  * fasta/fai/record: `Record::reference_sequence_name` is now `Record::name`.

    FASTA records are not necessarily reference sequences.

## 0.2.4 - 2021-09-23

### Fixed

  * fasta: Sync dependencies.

## 0.2.3 - 2021-09-19

### Fixed

  * fasta: Sync dependencies.

## 0.2.2 - 2021-08-19

### Fixed

  * fasta: Sync dependencies.

## 0.2.1 - 2021-08-11

### Fixed

  * fasta: Sync dependencies.

## 0.2.0 - 2021-07-30

### Changed

  * fasta/record: Rename `reference_sequence_name` to `name`.

    FASTA records are not necessarily reference sequences.

## 0.1.1 - 2021-07-21

### Fixed

  * fasta: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * fasta: Initial release.
