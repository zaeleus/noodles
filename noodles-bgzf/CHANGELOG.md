# Changelog

## Unreleased

### Changed

  * bgzf/writer: Immediately flush block upon reaching the max block size.

  * bgzf/writer: Track position (`Writer::position`) and virtual position
    (`Writer::virtual_position`).

## 0.9.0 - 2022-02-17

### Added

  * bgzf: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

  * bgzf/writer: Add builder (`Writer::builder`).

    This allows setting a compression level.

### Changed

  * bgzf: Update to tokio-util 0.7.0.

## 0.8.0 - 2022-01-27

### Changed

  * bgzf/async/reader: Create future for inflating data.

## 0.7.0 - 2021-12-02

### Added

  * bgzf/reader: Read to given buffer when it is guaranteed to be larger than
    the next block.

## 0.6.0 - 2021-11-18

### Added

  * bgzf/reader: Add common methods to access the underlying reader: `get_ref`,
    `get_mut`, and `into_inner`.

## 0.5.0 - 2021-11-11

### Added

  * bgzf: Add `libdeflate` feature to use [libdeflate] for encoding and
    decoding DEFLATE streams (#51).

  * bgzf/writer: Add crate `CompressionLevel` type.

    This replaces the usage of `flate2::Compression` and
    `libdeflater::CompressionLvl`.

[#51]: https://github.com/zaeleus/noodles/issues/51
[libdeflate]: https://github.com/ebiggers/libdeflate

### Changed

  * bgzf: Update to Rust 2021.

  * bgzf/async/writer/builder: The compression level wrapper changed from
    `flate2::Compression` to `noodles_bgzf::writer::CompressionLevel`.

## 0.4.0 - 2021-08-19

### Changed

  * bgzf: Update to tokio 1.10.0.

  * bgzf/async: I/O builders are now owned/consuming builders.

    This fixes the terminal method not being able to move out of a mutable
    reference.

### Fixed

  * bgzf: Define features to enable for Docs.rs.

## 0.3.0 - 2021-08-11

### Added

  * bgzf/async: Add async reader (`bgzf::AsyncReader`).

  * bgzf/async: Add async writer (`bgzf::AsyncWriter`) ([#17]).

    Async I/O can be enabled with the `async` feature. 

    Async BGZF I/O is implemented using a queue of block buffers, which are
    encoded/decoded in parallel (depending on the async executor). This can
    significantly improve read/write performance.

[#17]: https://github.com/zaeleus/noodles/issues/17

## 0.2.0 - 2021-07-21

### Fixed

  * bgzf: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * bgzf/index: Moved `Chunk` to noodles-csi.

    Replace usages of `noodles_bgzf::index::Chunk` with
    `noodles-csi::index::reference_sequence::bin::Chunk`.

  * bgzf/index: Moved `Metadata` to noodles-csi.

    Replace usages of `noodles_bgzf::index::Metadata` with
    `noodles-csi::index::reference_sequence::Metadata`.

  * bgzf/index: Moved chunk merging functions to noodles-csi.

    Replace usages of `noodles_bgzf::index::{merge_chunks, optimize_chunks}`
    with `noodles_csi::binning_index::{merge_chunks, optimize_chunks}`.

## 0.1.0 - 2021-07-14

  * bgzf: Initial release.
