# noodles

[![CI status](https://github.com/zaeleus/noodles/actions/workflows/ci.yml/badge.svg)](https://github.com/zaeleus/noodles/actions/workflows/ci.yml)

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers and writers for BAM 1.6, BCF 2.2, BGZF, CRAM 3.0,
FASTA, FASTQ, GFF3, SAM 1.6, tabix, and VCF 4.3.

Notably, the BAM and CRAM parsers are pure Rust implementations.

## Usage

There is currently no release of noodles, as the API is still experimental, but
it can be added to a test project as a preview.

noodles is split into multiple crates by file format. For convenience, a
top-level meta crate named `noodles` can be added to your project's dependency
list; and formats, listed as [features]. For example, to work with the BAM
format, enable the `bam` feature.

```toml
noodles = { git = "https://github.com/zaeleus/noodles.git", features = ["bam"] }
```

Each enabled feature can then be imported by its re-exported name, e.g.,

```rust
use noodles::bam;
```

[features]: https://doc.rust-lang.org/cargo/reference/features.html

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
