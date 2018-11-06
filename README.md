# noodles

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers (and some writers) for BAM, FASTQ, and GFF/GTFv2.

Notably, the BAM parser is a pure Rust implementation.

## Install

Add `noodles` to `Cargo.toml`.

```toml
noodles = { git = "https://github.com/zaeleus/noodles.git" }
```

## Related tools

noodles itself does not provide any applications, but related tools do depend
on it.

  * [noodles-count-features](https://github.com/zaeleus/noodles-count-features):
    Counts the number of aligned records that intersects a set of features.
