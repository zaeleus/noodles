# noodles

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers (and some writers) for BAM, FASTQ, and GFF/GTFv2.

Notably, the BAM parser is a pure Rust implementation.

## Related tools

noodles itself does not provide any applications, but related tools do depend
on it.

  * [noodles-count-features]: Counts the number of aligned records that
    intersects a set of features.

  * [noodles-fpkm]: Calculates FPKM values from feature counts.

[noodles-count-features]: https://github.com/zaeleus/noodles-count-features
[noodles-fpkm]: https://github.com/zaeleus/noodles-fpkm
