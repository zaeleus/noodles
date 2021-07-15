# Changelog

## Unreleased

### Added

  * index: Added shared `Metadata` struct.

    This is parsed from an optional pseudo-bin used in the binning index
    formats (BAI, CSI, and tabix). It contains start/end positions and
    mapped/unmapped record counts for a reference sequence.

## 0.1.0 - 2021-07-14

  * Initial release.
