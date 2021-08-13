# Changelog

## Unreleased

### Changed

  * Update to tokio 1.10.0.

### Fixed

  * Define features to enable for Docs.rs.

## 0.3.0 - 2021-08-11

### Added

  * async: Add async reader (`bgzf::AsyncReader`).

  * async: Add async writer (`bgzf::AsyncWriter`) ([#17]).

    Async I/O can be enabled with the `async` feature. 

    Async BGZF I/O is implemented using a queue of block buffers, which are
    encoded/decoded in parallel (depending on the async executor). This can
    signficantly improve read/write performance.

[#17]: https://github.com/zaeleus/noodles/issues/17

## 0.2.0 - 2021-07-21

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * index: Moved `Chunk` to noodles-csi.

    Replace usages of `noodles_bgzf::index::Chunk` with
    `noodles-csi::index::reference_sequence::bin::Chunk`.

  * index: Moved `Metadata` to noodles-csi.

    Replace usages of `noodles_bgzf::index::Metadata` with
    `noodles-csi::index::reference_sequence::Metadata`.

  * index: Moved chunk merging functions to noodles-csi.

    Replace usages of `noodles_bgzf::index::{merge_chunks, optimize_chunks}`
    with `noodles_csi::binning_index::{merge_chunks, optimize_chunks}`.

## 0.1.0 - 2021-07-14

  * Initial release.
