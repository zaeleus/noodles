# Changelog

## Unreleased

### Removed

  * bai/index/reference_sequence/metadata: Remove conversion to `Bin`.

    The metadata fields do not have a practical usage as bin chunks. The
    conversion was only used when writing the index, which now uses a metadata
    writer.

## 0.1.0 - 2021-07-14

  * Initial release.
