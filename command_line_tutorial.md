# noodles

[![GitHub Actions status](https://github.com/zaeleus/noodles/workflows/CI/badge.svg)](https://github.com/zaeleus/noodles/actions)

## Getting Started

In this tutorial, we will be creating a command line interface for Noodles similar to the [CLI Rust tutorial](https://rust-cli.github.io/book/index.html). At the end, you will have a working command line program which intakes the file, does something to it, and prints an output. 

Throughout Noodles, there are examples packaged with the source code for each file type (for example, noodles-bam has a src and an examples folder, which we will be using in this tutorial).

To get started, we first clone the noodles github in the folder we want to work in. For the purpose of this tutorial, we use the directory at `~/projects` and run:

`git clone https://github.com/zaeleus/noodles`

We should now have a folder in projects called noodles. In the noodles file, we have noodles-* and noodles. For this tutorial, we'll make a new folder called `noodles-tutorial`. We now should have `~/projects/noodles/noodles-tutorial`. This will be the folder we write our code in. We navigate to the `noodles-tutorial` folder and run `cargo init` to make a new cargo. At this point, we get a warning about adding our new folder as a member. We put the finishing touches on our workspace by adding the `noodles-tutorial` folder to the end of the members list in our top cargo.toml file, `~/projects/noodles/Cargo.toml`, which should now contain: 

`[workspace]`
`members = [`
`  "noodles",`
`  "noodles-bam",`
`  "noodles-bgzf",`
`  "noodles-cram",`
`  "noodles-fasta",`
`  "noodles-fastq",`
`  "noodles-gff",`
`  "noodles-sam",`
`  "noodles-tabix",`
`  "noodles-vcf",`
`  "noodles-tutorial",` `#we added it here`
`]`

We're all set to start working on our CLI. We make a new file in the `noodles-tutorial` directory and name it `main.rs`. It may be a little confusing, but another lower Cargo.toml file in `~/projects/noodles/noodles-tutorial/Cargo.toml` lists the dependencies we need for our cargo tutorial. We list all the noodles directories, along with the external crate for parsing our CLI arguments, `structopt`. The final cargo.toml in our noodles-tutorial now contains the dependencies:

`noodles = { path = "../noodles" }`
`noodles-bgzf = { path = "../noodles-bgzf" }`
`noodles-sam = { path = "../noodles-sam" }`
`noodles-bam = { path = "../noodles-bam" }`
`noodles-cram = { path = "../noodles-cram" }`
`noodles-fasta = { path = "../noodles-fasta" }`
`noodles-fastq = { path = "../noodles-fastq" }`
`noodles-gff = { path = "../noodles-gff" }`
`noodles-tabix = { path = "../noodles-tabix" }`
`structopt = "0.3.20"`

We fill our `main.rs` with one of the examples from noodles-sam, adding in our `infile` struct and parameter from the rust command line example.

We finally run `cargo build` in the `noodles-tutorial` directory, and we should see that the build finished (hopefully). We can now go to `~/projects/noodles/target/debug/`. There are a lot of files, but our CLI is the `noodles-tutorial` file. We can now run `./noodles-tutorial` and it will return that we didn't specify the input file. Cheers!