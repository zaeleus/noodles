# noodles

[![GitHub Actions status](https://github.com/zaeleus/noodles/workflows/CI/badge.svg)](https://github.com/zaeleus/noodles/actions)

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers (and some writers) for BAM, FASTQ, and GFF/GTFv2.

Notably, the BAM parser is a pure Rust implementation.

## Related tools

noodles itself does not provide any applications, but related tools do depend
on it.

  * [noodles-squab]: Gene expression quantification by counting aligned record
    intersections on reference gene annotations.

[noodles-squab]: https://github.com/zaeleus/noodles-squab
