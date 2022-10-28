# Changelog

## 0.17.0 - 2022-10-28

### Added

  * bgzf/async/reader/builder: Implement `Default`.

  * bgzf/async/writer/builder: Implement `Default`.

  * bgzf/writer/builder: Implement `Default`.

### Changed

  * bgzf/async/reader/builder: `Builder` no longer holds a reader.

  * bgzf/async/writer/builder: `Builder` no longer holds a reader.

  * bgzf/async/writer/builder: Change worker count to a `NonZeroUsize`.

    The worker count can never be 0.

  * bgzf/writer/builder: `Builder` no longer holds a writer.

### Removed

  * bgzf/async/reader: Remove `AsyncReader::builder`.

    Use `r#async::reader::Builder::default` instead.

  * bgzf/async/reader/builder: Remove `Builder::build`.

    Use `Builder::build_with_reader` instead.

  * bgzf/async/writer: Remove `AsyncWriter::builder`.

    Use `r#async::writer::Builder::default` instead.

  * bgzf/async/writer/builder: Remove `Builder::build`.

    Use `Builder::build_with_writer` instead.

  * bgzf/writer: Remove `Writer::builder`.

    Use `writer::Builder::default` instead.

  * bgzf/writer/builder: Remove `Builder::build`.

    Use `Builder::build_with_writer` instead.

## 0.16.0 - 2022-10-20

### Added

  * bgzf: Add an indexed reader (`IndexedReader`).

  * bgzf/reader/builder: Add build from path (`Builder::build_from_path`).

### Changed

  * bgzf: Update to libdeflater 0.11.0 (libdeflate 1.14).

## 0.15.0 - 2022-09-29

### Added

  * bgzf/gzi: Add convenience `read` function.

  * bgzf/reader: Add seeking to an uncompressed position
    (`Reader::seek_by_uncompressed_position`).

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.57.0.

  * bgzf/gzi: Start index with the initial block.

    The first block, i.e., (0, 0), is now included in the index.

### Fixed

  * bgzf/reader/block/multi: Reset block queue after possible inner reader
    mutation ([#112]).

[#112]: https://github.com/zaeleus/noodles/issues/112

## 0.14.0 - 2022-08-16

### Added

  * bgzf/gzi: Add async reader (`gzi::AsyncReader`).

  * bgzf/gzi: Add reader (`gzi::Reader`).

  * bgzf/reader: Add builder (`reader::Builder`).

  * bgzf/reader: Add multithreaded block reader.

    The multithreaded `bgzf::Reader` is enabled when the worker count is > 1.

### Changed

  * bgzf/async/reader/builder: Change worker count to a `NonZeroUsize`.

    The worker count can never be 0.

## 0.13.0 - 2022-07-05

### Changed

  * bgzf: Update to libdeflater 0.10.0 (libdeflate 1.10) ([#98]).

    This also fixes the usage of compression level 0 (e.g.,
    `CompressionLevel::none`) when the `libdeflate` feature is enabled.

  * bgzf/async/reader: Verify BGZF block header values ([#93]).

  * bgzf/reader: Verify BGZF block header values ([#93]).

[#93]: https://github.com/zaeleus/noodles/issues/93
[#98]: https://github.com/zaeleus/noodles/issues/98

### Fixed

  * bgzf/async/block_codec: Avoid truncating the header block size (`BSIZE`)
    ([#96]).

  * bgzf/async/writer: Compensate for gzip overhead ([#96]).

  * bgzf/writer: Avoid truncating the header block size (`BSIZE`) ([#96]).

  * bgzf/writer: Compensate for gzip overhead ([#96]).

    This reduces the max uncompressed data buffer to 65495 bytes and allows
    the use of compression level 0.

[#96]: https://github.com/zaeleus/noodles/issues/96

## 0.12.0 - 2022-06-08

### Changed

  * bgzf/async/reader: Verify the expected CRC32 checksum of the uncompressed
    data.

  * bgzf/reader: Verify the expected CRC32 checksum of the uncompressed data.

## 0.11.0 - 2022-03-29

### Added

  * bgzf/writer/compression_level: Enable compression level 0 (no compression)
    for libdeflate.

### Changed

  * bgzf: Update to libdeflater 0.8.0 (libdeflate 1.10).

## 0.10.0 - 2022-03-02

### Changed

  * bgzf/writer: Immediately flush block upon reaching the max block size.

  * bgzf/writer: Track position (`Writer::position`) and virtual position
    (`Writer::virtual_position`) ([#77]).

  * bgzf/writer: Add `into_inner` to return the underlying writer.

[#77]: https://github.com/zaeleus/noodles/issues/77

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
