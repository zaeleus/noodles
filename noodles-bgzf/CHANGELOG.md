# Changelog

## 0.42.0 - 2025-07-12

### Added

  * bgzf/async/fs: Add `open` function.

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.41.0 - 2025-05-29

### Removed

  * bgzf: Remove deprecated items.

    The following items are removed:

      * `AsyncReader` (deprecated since 0.35.0; use `r#async::io::Reader` instead),
      * `AsyncWriter` (0.35.0; `r#async::io::Writer`),
      * `IndexedReader` (0.38.0; `io::IndexedReader`),
      * `MultithreadedReader` (0.38.0; `io::MultithreadedReader`),
      * `MultithreadedWriter` (0.38.0; `io::MultithreadedWriter`),
      * `Reader` (0.38.0; `io::Reader`),
      * `Writer` (0.38.0; `io::Writer`),
      * `gzi::AsyncReader` (0.35.0; `gzi::r#async::io::Reader`),
      * `gzi::Reader` (0.35.0; `gzi::io::Reader`),
      * `gzi::r#async::Reader` (0.35.0; `gzi::r#async::io::Reader`),
      * `gzi::r#async::read` (0.35.0; `gzi::r#async::fs::read`),
      * `gzi::read` (0.35.0; `gzi::fs::read`),
      * `io::writer::Builder::build_with_writer` (0.33.0;
        `io::writer::Builder::build_from_writer`).
      * `io::writer::CompressionLevel::best` (0.29.0;
        `io::writer::CompressionLevel::BEST`),
      * `io::writer::CompressionLevel::fast` (0.29.0;
        `io::writer::CompressionLevel::FAST`),
      * `io::writer::CompressionLevel::none` (0.29.0;
        `io::writer::CompressionLevel::NONE`),
      * `r#async::Reader` (0.38.0; `r#async::io::Reader`),
      * `r#async::Writer` (0.38.0; `r#async::io::Writer`),
      * `r#async::io::reader::Builder::build_with_reader` (0.33.0;
        `r#async::io::reader::Builder::build_from_reader`), and
      * `r#async::io::writer::Builder::build_with_writer` (0.33.0;
        `r#async::io::writer::Builder::build_from_writer`),

## 0.40.0 - 2025-05-16

### Added

  * bgzf/fs: Add `open` function.

## 0.39.0 - 2025-04-13

### Fixed

  * bgzf/io/block/data: Remove buffer end mask ([#336]).

    The mask does not work when the buffer is equal to the max block size
    2**16.

[#336]: https://github.com/zaeleus/noodles/issues/336

## 0.38.0 - 2025-04-06

### Changed

  * bgzf: Move readers (`IndexedReader`, `MultithreadedReader`, and `Reader`)
    and writer (`MultithreadedWriter` and `Writer`) to `io` module.

### Deprecated

  * bgzf: Deprecate `IndexedReader`, `MultithreadedReader`,
    `MultithreadedWriter`, `Reader` and `Writer`.

    These moved to the `io` module.

## 0.37.0 - 2025-03-08

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.81.0.

  * bgzf: Change default DEFLATE backend to [zlib-rs] ([#331]).

[zlib-rs]: https://github.com/trifectatechfoundation/zlib-rs
[#331]: https://github.com/zaeleus/noodles/pull/331

## 0.36.0 - 2025-02-06

### Changed

  * bgzf: Sync dependencies.

## 0.35.0 - 2025-01-19

### Added

  * bgzf/gzi: Add wrapper for `Index`.

  * bgzf/gzi/fs: Add convenience write function (`write`).

  * bgzf/gzi/io: Add writer (`Writer`).

  * bgzf/gzi/async/io/reader: Add common methods to access the underlying I/O
    (`Reader::get_ref`, `Reader::get_mut`, and `Reader::into_inner`).

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.73.0.

  * bgzf/gzi: Move reader (`Reader`) to `io` module.

  * bgzf/gzi: Move convenience `read` function to `fs` module.

  * bgzf/gzi/index: Remove first entry.

    The first entry is now implicitly `(0, 0)`. As a result, the number of
    entries corresponds to n - 1 blocks. This now follows the same layout as
    the physical index.

### Deprecated

  * bgzf: Deprecate async re-exports (`AsyncReader` and `AsyncWriter`).

    Use `bgzf::r#async::Reader` and `bgzf::r#async::Writer` instead.

  * bgzf/gzi: Deprecate `Reader`.

    Use `bgzf::gzi::io::Reader` instead.

  * bgzf/gzi: Deprecate `AsyncReader` re-export.

    Use `bgzf::gzi::r#async::Reader` instead.

  * bgzf/gzi: Deprecate `read`.

    Use `gzi::fs::read` instead.

## 0.34.0 - 2024-12-12

### Added

  * bgzf/async/reader: Add position getter (`Reader::position`).

## 0.33.0 - 2024-09-04

### Changed

  * bgzf/reader/builder: Rename `Builder::build_with_reader` to
    `Builder::build_from_reader`.

  * bgzf/writer/builder: Rename `Builder::build_with_writer` to
    `Builder::build_from_writer`.

### Deprecated

  * bgzf/async/reader/builder: Deprecate `Builder::build_with_reader`.

    Use `Builder::build_from_reader` instead.

  * bgzf/writer/builder: Deprecate `Builder::build_with_writer`.

    Use `Builder::build_from_writer` instead.

## 0.32.0 - 2024-07-14

### Added

  * bgzf/async/reader: Add seeking to an uncompressed position
    (`Reader::seek_by_uncompressed_position`) ([#273]).

[#273]: https://github.com/zaeleus/noodles/issues/273

## 0.31.0 - 2024-06-17

### Added

  * bgzf/gzi/reader: Add common methods to access the underlying reader.

## 0.30.0 - 2024-05-16

### Added

  * bgzf/writer/compression_level: Add a `const` getter
    (`CompressionLevel::get`).

  * bgzf/writer/compression_level: Implement `Ord` + `PartialOrd`.

## 0.29.0 - 2024-05-02

### Added

  * bgzf/io: Add read-related traits (`Read`, `BufRead`, and `Seek`).

  * bgzf/multithreaded_reader: Add default constructor
    (`MultithreadedReader::new`) and accessor for the underlying reader
    (`MultithreadedReader::get_mut`).

  * bgzf/multithreaded_writer: Add default constructor
    (`MultithreadedWriter::new`).

  * bgzf/virtual_position: Add `const` constructor (`VirtualPosition::new`).

  * bgzf/writer/compression_level: Add associated constants for common
    compression levels (`CompressionLevel::NONE`, `CompressionLevel::FAST`, and
    `CompressionLevel::BEST`).

  * bgzf/writer/compression_level: Add `const` constructor
    (`CompressionLevel::new`).

### Changed

  * bgzf/multithreaded_writer: Return inner writer on shut down
    (`Multithreaded::finish`).

    This adds generic type parameter `W` to `MultithreadedWriter`.

  * bgzf/virtual_position: Make getters (`VirtualPosition::uncompressed` and
    `VirtualPosition::compressed`) `const`.

### Deprecated

  * bgzf/writer/compression_level: Deprecate named constructors
    (`CompressionLevel::none`, `CompressionLevel::fast`, and
    `CompressionLevel::best`).

    Use the compression level constants instead of the named constructors,
    e.g., `CompressionLevel::none()` => `CompressionLevel::NONE`.

### Removed

  * bgzf/reader/block: Remove blocking multithreaded reader.

    This block reader was used when the worker count was set to > 1. However,
    I/O was still performed on the current thread. In the case of sequential
    reading, use `bgzf::MultithreadedReader` instead, which moves I/O to its
    own thread.

## 0.28.0 - 2024-03-28

### Changed

  * bgzf/async/block_codec: Allow final frame to be incomplete ([#243]).

    The final frame of an input stream is now allowed to be incomplete,
    allowing for decoding up to a known end position in the stream. Frame size
    validation now occurs from the inflater.

[#243]: https://github.com/zaeleus/noodles/issues/243

## 0.27.0 - 2024-03-12

### Added

  * bgzf/multithreaded_writer: Add builder ([#238]).

    This also adds the ability to set the compression level.

[#238]: https://github.com/zaeleus/noodles/issues/238

## 0.26.0 - 2023-12-14

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.25.0 - 2023-10-12

### Changed

  * bgzf: Update to libdeflater 1.19.0 (libdeflate 1.19).

  * bgzf/multithreaded_reader: Return inner reader from
    `bgzf::MultithreadedReader::finish`.

## 0.24.0 - 2023-08-31

### Added

  * bgzf: Add a multithreaded reader (`bgzf::MultithreadedReader`).

    This differs from a `bgzf::Reader` with > 1 worker by placing the inner
    reader on its own thread to read raw frames asynchronously.

## 0.23.0 - 2023-08-17

### Added

  * bgzf/async/reader: Add inner access delegates (`r#async::Reader::get_ref`,
    `r#async::Reader::get_mut`, `r#async::Reader::get_pin_mut`, and
    `r#async::Reader::into_inner`).

## 0.22.0 - 2023-06-01

### Changed

  * bgzf: Update to libdeflater 0.14.0 (libdeflate 1.18).

## 0.21.0 - 2023-04-27

### Added

  * bgzf/indexed_reader: Add getter for index (`IndexedReader::index`).

### Changed

  * bgzf: Update to libdeflater 0.13.0.

## 0.20.0 - 2023-03-03

### Added

  * bgzf/gzi/async: Add convenience `read` function.

### Changed

  * bgzf/gzi: Improve performance of reading entire GZ index ([#143]).

[#143]: https://github.com/zaeleus/noodles/issues/143

## 0.19.0 - 2023-02-03

### Changed

  * bgzf: Raise minimum supported Rust version (MSRV) to 1.64.0.

  * bgzf: Update to libdeflater 0.12.0 (libdeflate 1.17).

### Removed

  * bgzf/virtual_position: Remove `VirtualPosition::max`.

    This was deprecated in noodles-bgzf 0.18.0. Use
    `VirtualPosition::MAX` instead.

## 0.18.0 - 2022-11-18

### Added

  * bgzf: Add multithreaded writer (`bgzf::MultithreadedWriter`).

    This is an alternative writer that uses a thread pool to compress block
    data. It currently cannot be used as a drop-in replacement for
    `bgzf::Writer`; however, it does implement `std::io::Write`.

  * bgzf/virtual_position: Add `MIN` and `MAX` associated constants.

### Fixed

  * bgzf/async/reader: Skip empty blocks.

  * bgzf/reader: Skip empty blocks ([#133]).

[#133]: https://github.com/zaeleus/noodles/issues/133

### Deprecated

  * bgzf/virtual_position: Deprecate `VirtualPosition::max`.

    Use `VirtualPosition::MAX` instead.

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
