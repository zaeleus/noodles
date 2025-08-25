# Changelog

## 0.57.0 - 2025-08-25

### Changed

  * tabix/io/writer/index/header/reference_sequence_names: Disallow `NUL` in
    names.

## 0.56.0 - 2025-07-12

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.55.0 - 2025-05-29

### Removed

  * tabix: Remove deprecated items.

    The following items are removed:

      * `AsyncReader` (deprecated since 0.45.0; use `r#async::io::Reader`
        instead),
      * `AsyncWriter` (0.45.0; `r#async::io::Writer`),
      * `Reader` (0.45.0; `io::Reader`),
      * `Writer` (0.45.0; `io::Writer`),
      * `r#async::Reader` (0.45.0; `r#async::io::Reader`),
      * `r#async::Writer` (0.45.0; `r#async::io::Writer`),
      * `r#async::read` (0.48.0; `r#async::fs::read`),
      * `r#async::write` (0.48.0; `r#async::fs::write`),
      * `read` (0.48.0; `fs::read`), and
      * `write` (0.48.0; `fs::write`).

## 0.54.0 - 2025-05-16

### Changed

  * tabix: Sync dependencies.

## 0.53.0 - 2025-04-13

### Changed

  * tabix: Sync dependencies.

## 0.52.0 - 2025-04-06

### Changed

  * tabix: Sync dependencies.

## 0.51.0 - 2025-03-08

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.50.0 - 2025-02-06

### Changed

  * tabix: Sync dependencies.

## 0.49.0 - 2025-01-23

### Changed

  * tabix: Sync dependencies.

## 0.48.0 - 2025-01-19

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.73.0.

  * tabix: Move convenience functions (`read` and `write`) to `fs` module.

### Deprecated

  * tabix: Deprecate `read` and `write`.

    Use `tabix::fs::read` and `tabix::fs::write`, respectively, instead.

## 0.47.0 - 2024-12-12

### Changed

  * tabix: Sync dependencies.

## 0.46.0 - 2024-11-07

### Changed

  * tabix: Sync dependencies.

## 0.45.0 - 2024-09-26

### Changed

  * tabix: Move reader (`Reader`) and writer (`Writer`) to `io` module.

  * tabix/async: Move reader (`Reader`) and writer (`Writer`) to `io`
    module.

### Deprecated

  * tabix: Deprecate `Reader` and `Writer`.

    Use `tabix::io::Reader` and `tabix::io::Writer`, respectively, instead.

  * tabix: Deprecate `AsyncReader` and `AsyncWriter`.

    Use `tabix::r#async::io::Reader` and `tabix::r#async::io::Writer`,
    respectively, instead.

## 0.44.0 - 2024-09-04

### Changed

  * tabix: Sync dependencies.

## 0.43.0 - 2024-07-14

### Changed

  * tabix: Update to bit-vec 0.7.0.

## 0.42.0 - 2024-06-17

### Added

  * tabix: Add common methods to access the underlying I/O.

### Changed

  * tabix/async/writer: `Writer::into_inner` now returns the inner BGZF
    writer instead of `R`.

    Use `writer.into_inner().into_inner()` to unwrap into `R`.

  * tabix/writer: `Writer::get_ref` now returns the inner BGZF writer
    instead of `R`.

    Use `writer.get_ref().get_ref()` to get `&R`.

## 0.41.0 - 2024-05-16

### Changed

  * tabix: Sync dependencies.

## 0.40.0 - 2024-05-08

### Changed

  * tabix: Sync dependencies.

## 0.39.0 - 2024-05-02

### Changed

  * tabix: Sync dependencies.

## 0.38.0 - 2024-03-28

### Changed

  * tabix: Sync dependencies.

## 0.37.0 - 2024-03-12

### Changed

  * tabix: Sync dependencies.

## 0.36.0 - 2024-01-25

### Changed

  * tabix: Sync dependencies.

## 0.35.0 - 2023-12-14

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * tabix: Define `tabix::Index` as `binning_index::Index<LinearIndex>`.

## 0.34.0 - 2023-11-14

### Added

  * tabix: Add `Index` type alias for `csi::Index`.

## 0.33.0 - 2023-11-13

### Changed

  * tabix: Sync dependencies.

## 0.32.0 - 2023-10-26

### Changed

  * tabix: Sync dependencies.

## 0.31.0 - 2023-10-19

### Changed

  * tabix: Sync dependencies.

## 0.30.0 - 2023-10-12

### Changed

  * tabix: Sync dependencies.

## 0.29.0 - 2023-08-31

### Changed

  * tabix: Sync dependencies.

## 0.28.0 - 2023-08-24

### Changed

  * tabix/async/reader: Disallow duplicate bin IDs.

## 0.27.0 - 2023-08-17

### Changed

  * tabix: Sync dependencies.

## 0.26.0 - 2023-08-03

### Changed

  * tabix/reader: Disallow duplicate bin IDs.

## 0.25.0 - 2023-07-06

### Changed

  * tabix: Sync dependencies.

## 0.24.0 - 2023-06-29

### Changed

  * tabix: Sync dependencies.

## 0.23.0 - 2023-06-15

### Changed

  * tabix: Sync dependencies.

## 0.22.0 - 2023-06-01

### Changed

  * tabix: Sync dependencies.

## 0.21.0 - 2023-05-18

### Added

  * tabix/io/indexed_reader: Add builder (`indexed_reader::Builder`).

### Changed

  * tabix/async/reader: Ensure header column indices are > 0.

  * tabix/reader: Ensure header column indices are > 0.

## 0.20.0 - 2023-05-04

### Changed

  * tabix: Sync dependencies.

## 0.19.0 - 2023-04-27

### Changed

  * tabix: Sync dependencies.

## 0.18.0 - 2023-04-06

### Changed

  * tabix: Replace `Index` with `noodles_csi::Index`.

### Removed

  * tabix/index: Remove `Header`.

    Use `noodles_csi::index::Header` instead.

## 0.17.0 - 2023-03-03

### Changed

  * tabix: Sync dependencies.

## 0.16.0 - 2023-02-03

### Added

  * tabix/index/header/format: Implement `std::error::Error::source` for
    `TryFromIntError`.

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * tabix/index: Remove `ReferenceSequenceNames`.

    This was deprecated in noodles-tabix 0.11.0. Use
    `noodles_tabix::index::header::ReferenceSequenceNames` instead.

  * tabix/index: Remove `Index::reference_sequence_names`.

    This was deprecated in noodles-tabix 0.11.0. Use
    `Header::reference_sequence_names` instead.

  * tabix/index: Remove `Index::unmapped_read_count`.

    This was deprecated in noodles-tabix 0.2.0. Use
    `Index::unplaced_unmapped_record_count` instead.

  * tabix/index/builder: Remove `Builder::set_unmapped_read_count`.

    This was deprecated in noodles-tabix 0.3.0. Use
    `Builder::set_unplaced_unmapped_record_count` instead.

## 0.15.0 - 2022-11-18

### Changed

  * tabix: Sync dependencies.

## 0.14.0 - 2022-10-28

### Changed

  * tabix: Sync dependencies.

## 0.13.0 - 2022-10-20

### Changed

  * tabix: Sync dependencies.

## 0.12.1 - 2022-09-29

### Fixed

  * tabix: Sync dependencies.

## 0.12.0 - 2022-08-16

### Changed

  * tabix: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.11.0 - 2022-07-05

### Changed

  * tabix/index/reference_sequence/bin: Change bin ID to a `usize`.

  * tabix/index: Move reference sequence names to `Header`.

### Deprecated

  * tabix/index: Deprecate `ReferenceSequenceNames`.

    Use `header::ReferenceSequenceNames` instead.

  * tabix/index: Deprecate `Index::reference_sequence_names`.

    Use `Header::reference_sequence_names` instead.

## 0.10.0 - 2022-06-09

### Changed

  * tabix/reader: Fail if reference sequence names buffer has trailing data.

### Fixed

  * tabix: Sync dependencies.

## 0.9.1 - 2022-06-08

### Fixed

  * tabix: Sync dependencies.

## 0.9.0 - 2022-03-29

### Changed

  * tabix/index/indexer: Change start and end to `Position`.

  * tabix/index/reference_sequence: `ReferenceSequence::query` returns an
    `io::Error` instead of `QueryError`.

### Fixed

  * tabix/index/reference_sequence: Ensure the start position is not out of
    range for a query (`2^29 - 1`).

## 0.8.1 - 2022-03-02

### Fixed

  * tabix: Sync dependencies.

## 0.8.0 - 2022-02-17

### Added

  * tabix: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.7.3 - 2022-01-27

### Fixed

  * tabix: Sync dependencies.

## 0.7.2 - 2021-12-02

### Fixed

  * tabix: Require tokio's `fs` feature as a dependency ([#62]).

[#62]: https://github.com/zaeleus/noodles/issues/62

## 0.7.1 - 2021-11-18

### Fixed

  * tabix: Sync dependencies.

## 0.7.0 - 2021-11-11

### Changed

  * tabix: Update to Rust 2021.

## 0.6.1 - 2021-09-19

### Fixed

  * tabix/async: Fix writer not finalizing.

## 0.6.0 - 2021-08-19

### Changed

  * tabix: Update to tokio 1.10.0.

### Fixed

  * tabix: Define features to enable for Docs.rs.

## 0.5.0 - 2021-08-11

### Added

  * tabix: Add convenience write function to write an index to a file:
    `tabix::write`.

  * tabix/async: Add async reader (`tabix::AsyncReader`).

  * tabix/async: Add async writer (`tabix::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.4.0 - 2021-08-04

### Changed

  * tabix/reader: Disallow duplicate reference sequence names.

## 0.3.0 - 2021-07-30

### Deprecated

  * tabix/index/builder: Deprecate `set_unmapped_read_count`.

    Use `set_unplaced_unmapped_record_count` instead.

### Fixed

  * tabix/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * tabix/index: Implemented `BinningIndex` for `Index`.

  * tabix/index: Added `query` method to find chunks that intersect the given
    region.

  * tabix/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

### Changed

  * tabix/index: Reference sequence names are stored as an
    `index::ReferenceSequenceNames` (`IndexSet<String>`).

### Deprecated

  * tabix/index: Deprecated `Index::unmapped_read_count`.

    Use `unplaced_unmapped_record_count` instead.

### Fixed

  * tabix: Fixed documentation link in package manifest ([#31]).

  * tabix/reader: Avoid casts that may truncate.

    Fields that convert from `i32` to other integer types now check whether
    they are in range.

  * tabix/writer: Avoid casts that may truncate.

    Fields that convert to `i32` from other integer types now check whether
    they are in range.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * tabix/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * tabix: Initial release.
