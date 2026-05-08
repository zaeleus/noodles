use std::io::{BufReader, Cursor};

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use noodles_fastq as fastq;

fn generate_fastq(num_records: usize, seq_len: usize) -> Vec<u8> {
    const BASES: &[u8] = b"ACGTACGTACGTACGT";
    const QUALS: &[u8] = b"FFFFFFFFFFFFFFF!";

    let record_size = 1 + 20 + 1 + seq_len + 1 + 2 + 1 + seq_len + 1;
    let mut data = Vec::with_capacity(num_records * record_size);

    for i in 0..num_records {
        data.extend_from_slice(b"@read_");
        data.extend_from_slice(format!("{i:012}").as_bytes());
        data.push(b'\n');

        for j in 0..seq_len {
            data.push(BASES[j % BASES.len()]);
        }
        data.push(b'\n');

        data.push(b'+');
        data.push(b'\n');

        for j in 0..seq_len {
            data.push(QUALS[j % QUALS.len()]);
        }
        data.push(b'\n');
    }

    data
}

fn bench_fastq(c: &mut Criterion, label: &str, num_records: usize, seq_len: usize) {
    let data = generate_fastq(num_records, seq_len);
    let data_len = data.len() as u64;

    let mut group = c.benchmark_group(format!("fastq/{label}"));
    group.throughput(Throughput::Bytes(data_len));

    // noodles: read_record (buffer reuse)
    group.bench_function(BenchmarkId::new("noodles_read_record", ""), |b| {
        b.iter(|| {
            let mut reader = fastq::io::Reader::new(BufReader::new(Cursor::new(&data)));
            let mut record = fastq::Record::default();
            let mut n = 0u64;
            let mut total_bases = 0u64;
            loop {
                match reader.read_record(&mut record) {
                    Ok(0) => break,
                    Ok(_) => {
                        n += 1;
                        total_bases += record.sequence().len() as u64;
                    }
                    Err(e) => panic!("{e}"),
                }
            }
            (n, total_bases)
        });
    });

    // noodles: records() iterator (clones per record)
    group.bench_function(BenchmarkId::new("noodles_records_iter", ""), |b| {
        b.iter(|| {
            let mut reader = fastq::io::Reader::new(BufReader::new(Cursor::new(&data)));
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

    // noodles: read_record via Builder (256k buffer, optimal for bulk scanning)
    group.bench_function(BenchmarkId::new("noodles_builder", ""), |b| {
        b.iter(|| {
            let mut reader = fastq::io::reader::Builder::default()
                .build_from_reader(BufReader::with_capacity(256 * 1024, Cursor::new(&data)));
            let mut record = fastq::Record::default();
            let mut n = 0u64;
            let mut total_bases = 0u64;
            loop {
                match reader.read_record(&mut record) {
                    Ok(0) => break,
                    Ok(_) => {
                        n += 1;
                        total_bases += record.sequence().len() as u64;
                    }
                    Err(e) => panic!("{e}"),
                }
            }
            (n, total_bases)
        });
    });

    // needletail
    group.bench_function(BenchmarkId::new("needletail", ""), |b| {
        b.iter(|| {
            let mut reader =
                needletail::parse_fastx_reader(Cursor::new(&data)).expect("valid FASTQ");
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
                helicase::FastqParser::<CONFIG, SliceInput<'_>>::from_slice(&data)
                    .expect("valid FASTQ");
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

fn bench_short_reads(c: &mut Criterion) {
    bench_fastq(c, "short_reads_150bp_100k", 100_000, 150);
}

fn bench_long_reads(c: &mut Criterion) {
    bench_fastq(c, "long_reads_10kbp_10k", 10_000, 10_000);
}

criterion_group!(benches, bench_short_reads, bench_long_reads);
criterion_main!(benches);
