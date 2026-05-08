use std::io::{BufReader, Cursor};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use noodles_fasta as fasta;

fn generate_fasta(num_records: usize, seq_len: usize, line_width: usize) -> Vec<u8> {
    const BASES: &[u8] = b"ACGTACGTACGTACGT";

    let num_lines = (seq_len + line_width - 1) / line_width;
    let record_size = 1 + 20 + 1 + seq_len + num_lines + 1;
    let mut data = Vec::with_capacity(num_records * record_size);

    for i in 0..num_records {
        data.extend_from_slice(b">seq_");
        data.extend_from_slice(format!("{i:012}").as_bytes());
        data.push(b'\n');

        for j in 0..seq_len {
            if j > 0 && j % line_width == 0 {
                data.push(b'\n');
            }
            data.push(BASES[j % BASES.len()]);
        }
        data.push(b'\n');
    }

    data
}

fn bench_fasta(
    c: &mut Criterion,
    label: &str,
    num_records: usize,
    seq_len: usize,
    line_width: usize,
) {
    let data = generate_fasta(num_records, seq_len, line_width);
    let data_len = data.len() as u64;

    let mut group = c.benchmark_group(format!("fasta/{label}"));
    group.throughput(Throughput::Bytes(data_len));

    // noodles: records() iterator
    group.bench_function(BenchmarkId::new("noodles_records_iter", ""), |b| {
        b.iter(|| {
            let mut reader = fasta::io::Reader::new(BufReader::new(Cursor::new(&data)));
            let mut n = 0u64;
            let mut total_bases = 0u64;
            for result in reader.records() {
                let record = result.unwrap();
                n += 1;
                total_bases += record.sequence().len() as u64;
            }
            (n, total_bases)
        });
    });

    // needletail
    group.bench_function(BenchmarkId::new("needletail", ""), |b| {
        b.iter(|| {
            let mut reader =
                needletail::parse_fastx_reader(Cursor::new(&data)).expect("valid FASTA");
            let mut n = 0u64;
            let mut total_bases = 0u64;
            while let Some(result) = reader.next() {
                let record = result.unwrap();
                n += 1;
                total_bases += record.seq().len() as u64;
            }
            (n, total_bases)
        });
    });

    // helicase (from slice — in-memory, SIMD-accelerated)
    group.bench_function(BenchmarkId::new("helicase", ""), |b| {
        use helicase::input::{FromSlice, SliceInput};
        use helicase::parser::Event;
        use helicase::{Config, HelicaseParser, ParserOptions};

        const CONFIG: Config = ParserOptions::default().config();

        b.iter(|| {
            let mut parser =
                helicase::FastaParser::<CONFIG, SliceInput<'_>>::from_slice(&data)
                    .expect("valid FASTA");
            let mut n = 0u64;
            let mut total_bases = 0u64;
            while let Some(event) = parser.next() {
                match event {
                    Event::Record(_) => {
                        n += 1;
                        total_bases += parser.get_dna_string().len() as u64;
                    }
                    _ => {}
                }
            }
            (n, total_bases)
        });
    });

    group.finish();
}

fn bench_many_short_seqs(c: &mut Criterion) {
    bench_fasta(c, "short_1kbp_10k_w80", 10_000, 1_000, 80);
}

fn bench_few_long_seqs(c: &mut Criterion) {
    bench_fasta(c, "long_1mbp_25_w80", 25, 1_000_000, 80);
}

fn bench_single_line(c: &mut Criterion) {
    bench_fasta(c, "short_1kbp_10k_single_line", 10_000, 1_000, usize::MAX);
}

criterion_group!(benches, bench_many_short_seqs, bench_few_long_seqs, bench_single_line);
criterion_main!(benches);
