# noodles fastx

A FASTA and FASTQ parser.

## Evaluate parser

### Criterion

Requirement:
- [cargo criterion](https://github.com/bheisler/cargo-criterion)

```
cargo criterion
```

A html report is generate in `../target/criterion/reports/index.html`

### Hyperfine

Requirement:
- [hyperfine](https://github.com/sharkdp/hyperfine/)
- [seqtk](https://github.com/lh3/seqtk)
- a fastq file 

Next script assume variabale FASTQ contain path to fastq file.

```
cargo build --release --example fastq2fasta
hyperfine --warmup 3 -n fastx -n seqtk '../target/release/examples/fastq2fasta $FASTQ ' 'seqtk seq -A $FASTQ '
```
