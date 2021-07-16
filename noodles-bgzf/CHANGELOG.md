# Changelog

## Unreleased

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
