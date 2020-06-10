# noodles

[![GitHub Actions status](https://github.com/zaeleus/noodles/workflows/CI/badge.svg)](https://github.com/zaeleus/noodles/actions)

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers (and some writers) for BAM, BGZF, CRAM 3.0, FASTA,
FASTQ, GFF/GTFv2, SAM, and VCF 4.3.

Notably, the BAM and CRAM parsers are pure Rust implementations.

## Usage

There is currently no release of noodles, as the API is still experimental, but
it can be added to a test project as a preview.

noodles is split into multiple crates by file format. For example, to work with
the BAM format, add `noodles-bam` to your project's dependencies list.

```toml
noodles-bam = { git = "https://github.com/zaeleus/noodles.git" }
```

## Examples

Each crate may have its own examples directory, and all examples are runnable
as an application. After cloning the repository, run `cargo run --release
--example` for a list of available examples. Use the example name as the option
argument and append program arguments to the command, e.g.,

```bash
cargo run --release --example bam_write > sample.bam
cargo run --release --example bam_read_header sample.bam
```

## Related tools

noodles itself does not provide any applications, but related tools do depend
on it.

  * [noodles-squab]: Gene expression quantification by counting aligned record
    intersections on reference gene annotations.

[noodles-squab]: https://github.com/zaeleus/noodles-squab
