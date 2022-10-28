# noodles

[![crates.io](https://img.shields.io/crates/v/noodles.svg)](https://crates.io/crates/noodles)
[![Docs.rs](https://docs.rs/noodles/badge.svg)](https://docs.rs/noodles)
[![CI status](https://github.com/zaeleus/noodles/actions/workflows/ci.yml/badge.svg)](https://github.com/zaeleus/noodles/actions/workflows/ci.yml)

**noodles** is a library for handling various bioinformatics file formats. It
currently includes readers and writers for BAM 1.6, BCF 2.2, BED, BGZF, CRAM
3.0, CSI, FASTA, FASTQ, GFF3, GTF 2.2, SAM 1.6, tabix, and VCF 4.3.

## Usage

noodles is published on [crates.io]. Early versions can be used in projects,
but keep in mind that the API is still considered experimental.

noodles is split into multiple crates by file format. For convenience, a
top-level meta crate named `noodles` can be added to your project's dependency
list; and formats, listed as [features]. For example, to work with the BAM
format, enable the `bam` feature.

```toml
noodles = { version = "0.29.0", features = ["bam"] }
```

Each enabled feature can then be imported by its re-exported name, e.g.,

```rust
use noodles::bam;
```

[crates.io]: https://crates.io/
[features]: https://doc.rust-lang.org/cargo/reference/features.html

### Feature flags

Individual crates may have optional features that can be enabled using feature
flags.

  * `async`: Enables asynchronous I/O with [Tokio]. (BAM, BCF, BGZF, CRAM, CSI,
    FASTA, FASTQ, SAM, tabix, and VCF)
  * `libdeflate`: Use [libdeflate] to encode and decode DEFLATE streams. (BGZF)

[Tokio]: https://tokio.rs/
[libdeflate]: https://github.com/ebiggers/libdeflate

## Examples

Each crate may have its own examples directory, and all examples are runnable
as an application. After cloning the repository, run `cargo run --release
--example` for a list of available examples. Use the example name as the option
argument and append program arguments to the command, e.g.,

```bash
cargo run --release --example bam_write > sample.bam
cargo run --release --example bam_read_header sample.bam
```
