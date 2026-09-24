//! Benchmarks for IndexedReader::query performance.
//!
//! Run with: `cargo bench -p noodles-fasta`
//!
//! For the hs38DH-like benchmark:
//!   cargo bench -p noodles-fasta --bench indexed_reader -- "hs38dh"

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Cursor, Write};
use std::path::PathBuf;

use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};

use noodles_core::{Position, Region};
use noodles_fasta::{self as fasta, fai};

/// Path to the hs38DH FAI file for structure reference.
const HS38DH_FAI_PATH: &str = "/Users/nhomer/work/references/hs38DH/hs38DH.fa.fai";

/// Generate a FASTA sequence with the given length and line width.
fn generate_fasta(name: &str, seq_len: usize, line_bases: usize) -> Vec<u8> {
    let mut fasta = Vec::new();

    // Header
    fasta.push(b'>');
    fasta.extend_from_slice(name.as_bytes());
    fasta.push(b'\n');

    // Sequence with line breaks
    let bases = b"ACGTACGTACGTACGT";
    for i in 0..seq_len {
        fasta.push(bases[i % bases.len()]);
        if (i + 1) % line_bases == 0 {
            fasta.push(b'\n');
        }
    }
    // Final newline if not already added
    if seq_len % line_bases != 0 {
        fasta.push(b'\n');
    }

    fasta
}

/// Create an FAI record for the generated FASTA.
fn create_fai_record(name: &str, seq_len: usize, line_bases: usize) -> fai::Record {
    let header_len = 1 + name.len() + 1; // '>' + name + '\n'
    let line_width = line_bases + 1; // bases + '\n'
    fai::Record::new(
        name,
        seq_len as u64,
        header_len as u64,
        line_bases as u64,
        line_width as u64,
    )
}

/// Parse FAI records from a file (name, length, line_bases, line_width).
fn parse_fai_file(path: &str) -> Vec<(String, usize, usize, usize)> {
    let file = File::open(path).expect("Failed to open FAI file");
    let reader = BufReader::new(file);

    reader
        .lines()
        .map(|line| {
            let line = line.expect("Failed to read line");
            let fields: Vec<&str> = line.split('\t').collect();
            let name = fields[0].to_string();
            let length: usize = fields[1].parse().expect("Invalid length");
            let line_bases: usize = fields[3].parse().expect("Invalid line_bases");
            let line_width: usize = fields[4].parse().expect("Invalid line_width");
            (name, length, line_bases, line_width)
        })
        .collect()
}

/// Generate a synthetic FASTA file matching hs38DH structure.
/// Returns the path to the generated file and the FAI index.
fn generate_hs38dh_like_fasta() -> (PathBuf, fai::Index) {
    let fai_records = parse_fai_file(HS38DH_FAI_PATH);

    // Use a temp file
    let fasta_path = std::env::temp_dir().join("bench_hs38dh_like.fa");

    eprintln!(
        "Generating synthetic FASTA with {} contigs at {:?}...",
        fai_records.len(),
        fasta_path
    );

    let file = File::create(&fasta_path).expect("Failed to create temp FASTA");
    let mut writer = BufWriter::with_capacity(64 * 1024 * 1024, file);

    let mut index_records = Vec::with_capacity(fai_records.len());
    let mut offset: u64 = 0;

    let bases = b"ACGTACGTACGTACGT";

    for (name, length, line_bases, line_width) in &fai_records {
        // Write header
        write!(writer, ">{}\n", name).unwrap();
        let header_len = 1 + name.len() + 1;
        let seq_offset = offset + header_len as u64;

        // Write sequence with proper line breaks
        let mut pos = 0;
        while pos < *length {
            let end = (*length).min(pos + line_bases);
            for i in pos..end {
                writer.write_all(&[bases[i % bases.len()]]).unwrap();
            }
            writer.write_all(b"\n").unwrap();
            pos = end;
        }

        // Calculate total bytes for this record
        let num_full_lines = length / line_bases;
        let last_line_bases = length % line_bases;
        let seq_bytes = if last_line_bases > 0 {
            num_full_lines * line_width + last_line_bases + 1
        } else {
            num_full_lines * line_width
        };

        index_records.push(fai::Record::new(
            name.as_str(),
            *length as u64,
            seq_offset,
            *line_bases as u64,
            *line_width as u64,
        ));

        offset += header_len as u64 + seq_bytes as u64;
    }

    writer.flush().unwrap();
    drop(writer);

    let file_size = std::fs::metadata(&fasta_path).unwrap().len();
    eprintln!(
        "Generated {:.2} GB FASTA file",
        file_size as f64 / (1024.0 * 1024.0 * 1024.0)
    );

    (fasta_path, fai::Index::from(index_records))
}

/// Benchmark with hs38DH-like structure (real file I/O).
fn bench_hs38dh_like(c: &mut Criterion) {
    // Check if FAI file exists
    if !std::path::Path::new(HS38DH_FAI_PATH).exists() {
        eprintln!(
            "Skipping hs38dh benchmark: FAI file not found at {}",
            HS38DH_FAI_PATH
        );
        return;
    }

    let mut group = c.benchmark_group("indexed_reader_hs38dh_like");
    group.sample_size(20); // Fewer samples for slower benchmarks

    let (fasta_path, index) = generate_hs38dh_like_fasta();

    // Get contig info for queries
    let fai_records = parse_fai_file(HS38DH_FAI_PATH);
    let large_contigs: Vec<_> = fai_records
        .iter()
        .filter(|(_, len, _, _)| *len > 10_000_000) // Only contigs > 10MB
        .collect();

    // Benchmark: Many small queries across large contigs
    let query_size = 1000;
    let num_queries = 10_000;

    group.throughput(Throughput::Elements(num_queries as u64));

    group.bench_function(
        BenchmarkId::new("random_queries", format!("{num_queries}x{query_size}bp")),
        |b| {
            b.iter(|| {
                let file = File::open(&fasta_path).unwrap();
                let buf_reader = BufReader::new(file);
                let mut reader = fasta::io::IndexedReader::new(buf_reader, index.clone());

                for i in 0..num_queries {
                    // Pick a random large contig
                    let (name, length, _, _) = large_contigs[i % large_contigs.len()];
                    let max_start = length - query_size;
                    let start = (i * 997) % max_start;

                    let start_pos = Position::try_from(start + 1).unwrap();
                    let end_pos = Position::try_from(start + query_size).unwrap();
                    let region = Region::new(name.as_str(), start_pos..=end_pos);
                    let record = reader.query(&region).unwrap();
                    black_box(record);
                }
            });
        },
    );

    // Benchmark: Larger queries
    let query_size = 100_000;
    let num_queries = 1_000;

    group.throughput(Throughput::Bytes((num_queries * query_size) as u64));

    group.bench_function(
        BenchmarkId::new("large_queries", format!("{num_queries}x{query_size}bp")),
        |b| {
            b.iter(|| {
                let file = File::open(&fasta_path).unwrap();
                let buf_reader = BufReader::new(file);
                let mut reader = fasta::io::IndexedReader::new(buf_reader, index.clone());

                for i in 0..num_queries {
                    let (name, length, _, _) = large_contigs[i % large_contigs.len()];
                    let max_start = length - query_size;
                    let start = (i * 997) % max_start;

                    let start_pos = Position::try_from(start + 1).unwrap();
                    let end_pos = Position::try_from(start + query_size).unwrap();
                    let region = Region::new(name.as_str(), start_pos..=end_pos);
                    let record = reader.query(&region).unwrap();
                    black_box(record);
                }
            });
        },
    );

    // Benchmark: Single-base queries (like fgumi's review command)
    let query_size = 1;
    let num_queries = 100_000;

    group.throughput(Throughput::Elements(num_queries as u64));

    group.bench_function(
        BenchmarkId::new(
            "single_base_queries",
            format!("{num_queries}x{query_size}bp"),
        ),
        |b| {
            b.iter(|| {
                let file = File::open(&fasta_path).unwrap();
                let buf_reader = BufReader::new(file);
                let mut reader = fasta::io::IndexedReader::new(buf_reader, index.clone());

                for i in 0..num_queries {
                    let (name, length, _, _) = large_contigs[i % large_contigs.len()];
                    let pos = (i * 997) % (length - 1);

                    let start_pos = Position::try_from(pos + 1).unwrap();
                    let region = Region::new(name.as_str(), start_pos..=start_pos);
                    let record = reader.query(&region).unwrap();
                    black_box(record);
                }
            });
        },
    );

    // Benchmark: Load ALL contigs using LINE-BY-LINE reader (original fgumi approach)
    // This was the original approach before using FAI index
    group.bench_function(
        BenchmarkId::new("load_all_contigs_line_by_line", "3366_contigs"),
        |b| {
            b.iter(|| {
                let file = File::open(&fasta_path).unwrap();
                let buf_reader = BufReader::new(file);
                let mut reader = fasta::io::Reader::new(buf_reader);

                let mut total_bases = 0usize;
                for result in reader.records() {
                    let record = result.unwrap();
                    total_bases += record.sequence().len();
                    black_box(&record);
                }
                black_box(total_bases);
            });
        },
    );

    // Benchmark: Load ALL contigs using IndexedReader (optimized approach)
    group.bench_function(
        BenchmarkId::new("load_all_contigs_indexed", "3366_contigs"),
        |b| {
            b.iter(|| {
                let file = File::open(&fasta_path).unwrap();
                let buf_reader = BufReader::new(file);
                let mut reader = fasta::io::IndexedReader::new(buf_reader, index.clone());

                let mut total_bases = 0usize;
                for (name, length, _, _) in &fai_records {
                    let start_pos = Position::try_from(1_usize).unwrap();
                    let end_pos = Position::try_from(*length).unwrap();
                    let region = Region::new(name.as_str(), start_pos..=end_pos);
                    let record = reader.query(&region).unwrap();
                    total_bases += record.sequence().len();
                    black_box(&record);
                }
                black_box(total_bases);
            });
        },
    );

    group.finish();

    // Cleanup
    std::fs::remove_file(&fasta_path).ok();
}

/// Benchmark many small queries (tests per-query overhead).
fn bench_many_small_queries(c: &mut Criterion) {
    let mut group = c.benchmark_group("indexed_reader_small_queries");

    let seq_len = 1_000_000; // 1 MB sequence
    let line_bases = 70; // Match hs38DH
    let query_size = 100; // 100 bp per query
    let num_queries = 10_000;

    let fasta_data = generate_fasta("chr1", seq_len, line_bases);
    let fai_record = create_fai_record("chr1", seq_len, line_bases);
    let index = fai::Index::from(vec![fai_record]);

    group.throughput(Throughput::Elements(num_queries as u64));

    group.bench_function(
        BenchmarkId::new("query", format!("{num_queries}x{query_size}bp")),
        |b| {
            b.iter(|| {
                let mut reader =
                    fasta::io::IndexedReader::new(Cursor::new(&fasta_data), index.clone());

                for i in 0..num_queries {
                    let start = (i * 97) % (seq_len - query_size); // Pseudo-random positions
                    let start_pos = Position::try_from(start + 1).unwrap();
                    let end_pos = Position::try_from(start + query_size).unwrap();
                    let region = Region::new("chr1", start_pos..=end_pos);
                    let record = reader.query(&region).unwrap();
                    black_box(record);
                }
            });
        },
    );

    group.finish();
}

/// Benchmark queries of varying sizes (tests scaling with query size).
fn bench_query_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("indexed_reader_query_sizes");

    let seq_len = 10_000_000; // 10 MB sequence
    let line_bases = 70; // Match hs38DH

    let fasta_data = generate_fasta("chr1", seq_len, line_bases);
    let fai_record = create_fai_record("chr1", seq_len, line_bases);
    let index = fai::Index::from(vec![fai_record]);

    // Test different query sizes
    for query_size in [100, 1_000, 10_000, 100_000, 1_000_000] {
        let num_queries = 100;

        group.throughput(Throughput::Bytes((query_size * num_queries) as u64));

        group.bench_with_input(
            BenchmarkId::new("query", format!("{query_size}bp")),
            &query_size,
            |b, &query_size| {
                b.iter(|| {
                    let mut reader =
                        fasta::io::IndexedReader::new(Cursor::new(&fasta_data), index.clone());

                    for i in 0..num_queries {
                        let start = (i * 1000) % (seq_len - query_size);
                        let start_pos = Position::try_from(start + 1).unwrap();
                        let end_pos = Position::try_from(start + query_size).unwrap();
                        let region = Region::new("chr1", start_pos..=end_pos);
                        let record = reader.query(&region).unwrap();
                        black_box(record);
                    }
                });
            },
        );
    }

    group.finish();
}

/// Benchmark single-line queries (tests fast path).
fn bench_single_line_queries(c: &mut Criterion) {
    let mut group = c.benchmark_group("indexed_reader_single_line");

    let seq_len = 1_000_000;
    let line_bases = 70; // Match hs38DH
    let num_queries = 10_000;

    let fasta_data = generate_fasta("chr1", seq_len, line_bases);
    let fai_record = create_fai_record("chr1", seq_len, line_bases);
    let index = fai::Index::from(vec![fai_record]);

    // Query sizes that fit within a single line (70bp lines)
    // Include 1bp to match fgumi's review command pattern
    for query_size in [1, 10, 35, 69] {
        group.throughput(Throughput::Elements(num_queries as u64));

        group.bench_with_input(
            BenchmarkId::new("query", format!("{query_size}bp")),
            &query_size,
            |b, &query_size| {
                b.iter(|| {
                    let mut reader =
                        fasta::io::IndexedReader::new(Cursor::new(&fasta_data), index.clone());

                    for i in 0..num_queries {
                        // Ensure we stay within a single line
                        let line_num = (i * 7) % (seq_len / line_bases);
                        let start = line_num * line_bases;
                        let start_pos = Position::try_from(start + 1).unwrap();
                        let end_pos = Position::try_from(start + query_size).unwrap();
                        let region = Region::new("chr1", start_pos..=end_pos);
                        let record = reader.query(&region).unwrap();
                        black_box(record);
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_hs38dh_like,
    bench_many_small_queries,
    bench_query_sizes,
    bench_single_line_queries,
);
criterion_main!(benches);
