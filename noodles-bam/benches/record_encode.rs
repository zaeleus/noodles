//! Benchmarks for BAM record encoding.
//!
//! This module benchmarks the performance of various BAM record encoding strategies,
//! comparing baseline implementations against optimized versions.

use std::io;

use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};

use noodles_bam::{
    self as bam,
    record::codec::encoder::{encode_record_buf, encode_with_prealloc, estimate_record_size},
};
use noodles_sam::{
    self as sam,
    alignment::{
        RecordBuf,
        record::{
            Flags, MappingQuality,
            cigar::{Op, op::Kind},
        },
        record_buf::{QualityScores, Sequence},
    },
};

/// Creates a realistic test record with the specified sequence length.
fn create_test_record(seq_len: usize) -> RecordBuf {
    let sequence: Vec<u8> = (0..seq_len)
        .map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        })
        .collect();

    let quality_scores: Vec<u8> = (0..seq_len).map(|i| ((i % 42) + 10) as u8).collect();

    RecordBuf::builder()
        .set_name("test_read_with_longer_name")
        .set_flags(Flags::SEGMENTED | Flags::FIRST_SEGMENT)
        .set_mapping_quality(MappingQuality::new(30).expect("valid mapping quality"))
        .set_cigar([Op::new(Kind::Match, seq_len)].into_iter().collect())
        .set_sequence(Sequence::from(sequence))
        .set_quality_scores(QualityScores::from(quality_scores))
        .build()
}

fn bench_encode_record(c: &mut Criterion) {
    let header = sam::Header::default();

    let mut group = c.benchmark_group("encode_record");

    // Test various sequence lengths
    for seq_len in [100, 150, 300, 500, 1000] {
        let record = create_test_record(seq_len);

        group.throughput(Throughput::Elements(1));

        // Benchmark with fresh buffer each time (simulates worst case)
        group.bench_with_input(
            BenchmarkId::new("fresh_buffer", seq_len),
            &seq_len,
            |b, _| {
                b.iter(|| {
                    let mut buf = Vec::new();
                    encode_with_prealloc(&mut buf, &header, black_box(&record)).unwrap();
                    black_box(buf)
                });
            },
        );

        // Benchmark with reused buffer (simulates typical batch encoding)
        group.bench_with_input(
            BenchmarkId::new("reused_buffer", seq_len),
            &seq_len,
            |b, _| {
                let mut buf = Vec::with_capacity(1024);
                b.iter(|| {
                    buf.clear();
                    encode_with_prealloc(&mut buf, &header, black_box(&record)).unwrap();
                    buf.len()
                });
                black_box(&buf);
            },
        );
    }

    group.finish();
}

fn bench_estimate_record_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("estimate_record_size");

    for seq_len in [100, 150, 300, 500, 1000] {
        let record = create_test_record(seq_len);

        group.bench_with_input(BenchmarkId::new("estimate", seq_len), &seq_len, |b, _| {
            b.iter(|| black_box(estimate_record_size(black_box(&record))));
        });
    }

    group.finish();
}

fn bench_batch_encode(c: &mut Criterion) {
    let header = sam::Header::default();

    let mut group = c.benchmark_group("batch_encode");

    // Create a batch of records
    let batch_size = 1000;
    let records: Vec<RecordBuf> = (0..batch_size)
        .map(|i| create_test_record(100 + (i % 100)))
        .collect();

    let total_seq_len: usize = records.iter().map(|r| r.sequence().len()).sum();

    group.throughput(Throughput::Bytes(total_seq_len as u64));

    // Benchmark batch encoding with pre-allocation
    group.bench_function("with_prealloc", |b| {
        let mut buf = Vec::new();
        b.iter(|| {
            buf.clear();
            for record in &records {
                encode_with_prealloc(&mut buf, &header, black_box(record)).unwrap();
            }
            buf.len()
        });
        black_box(&buf);
    });

    group.finish();
}

fn bench_encode_record_buf_vs_generic(c: &mut Criterion) {
    let header = sam::Header::default();

    let mut group = c.benchmark_group("encode_comparison");

    // Test typical read lengths
    for seq_len in [100, 150, 300] {
        let record = create_test_record(seq_len);

        group.throughput(Throughput::Elements(1));

        // Generic encoder with prealloc
        group.bench_with_input(
            BenchmarkId::new("generic_prealloc", seq_len),
            &seq_len,
            |b, _| {
                let mut buf = Vec::with_capacity(1024);
                b.iter(|| {
                    buf.clear();
                    encode_with_prealloc(&mut buf, &header, black_box(&record)).unwrap();
                    buf.len()
                });
                black_box(&buf);
            },
        );

        // Optimized RecordBuf encoder
        group.bench_with_input(
            BenchmarkId::new("record_buf_optimized", seq_len),
            &seq_len,
            |b, _| {
                let mut buf = Vec::with_capacity(1024);
                b.iter(|| {
                    buf.clear();
                    encode_record_buf(&mut buf, &header, black_box(&record)).unwrap();
                    buf.len()
                });
                black_box(&buf);
            },
        );
    }

    group.finish();
}

fn bench_batch_encode_comparison(c: &mut Criterion) {
    let header = sam::Header::default();

    let mut group = c.benchmark_group("batch_encode_comparison");

    // Create a batch of records
    let batch_size = 1000;
    let records: Vec<RecordBuf> = (0..batch_size)
        .map(|i| create_test_record(100 + (i % 100)))
        .collect();

    let total_seq_len: usize = records.iter().map(|r| r.sequence().len()).sum();

    group.throughput(Throughput::Bytes(total_seq_len as u64));

    // Batch encoding with generic encoder
    group.bench_function("generic_prealloc", |b| {
        let mut buf = Vec::new();
        b.iter(|| {
            buf.clear();
            for record in &records {
                encode_with_prealloc(&mut buf, &header, black_box(record)).unwrap();
            }
            buf.len()
        });
        black_box(&buf);
    });

    // Batch encoding with optimized encoder
    group.bench_function("record_buf_optimized", |b| {
        let mut buf = Vec::new();
        b.iter(|| {
            buf.clear();
            for record in &records {
                encode_record_buf(&mut buf, &header, black_box(record)).unwrap();
            }
            buf.len()
        });
        black_box(&buf);
    });

    group.finish();
}

fn bench_writer_methods(c: &mut Criterion) {
    use sam::alignment::io::Write as _;

    let header = sam::Header::default();

    let mut group = c.benchmark_group("writer_methods");

    // Create a batch of records
    let batch_size = 1000;
    let records: Vec<RecordBuf> = (0..batch_size)
        .map(|i| create_test_record(100 + (i % 100)))
        .collect();

    let total_seq_len: usize = records.iter().map(|r| r.sequence().len()).sum();

    group.throughput(Throughput::Bytes(total_seq_len as u64));

    // Writer with generic write_alignment_record
    group.bench_function("write_alignment_record", |b| {
        b.iter(|| {
            let mut writer = bam::io::Writer::from(io::sink());
            for record in &records {
                writer
                    .write_alignment_record(&header, black_box(record))
                    .unwrap();
            }
        });
    });

    // Writer with optimized write_record_buf
    group.bench_function("write_record_buf", |b| {
        b.iter(|| {
            let mut writer = bam::io::Writer::from(io::sink());
            for record in &records {
                writer.write_record_buf(&header, black_box(record)).unwrap();
            }
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_encode_record,
    bench_estimate_record_size,
    bench_batch_encode,
    bench_encode_record_buf_vs_generic,
    bench_batch_encode_comparison,
    bench_writer_methods,
);

criterion_main!(benches);
