/* benchmark use */
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::seq::SliceRandom;

/* worker use */
use noodles_fasta;
use noodles_fastx;

fn generate_fasta(nb_reads: usize, nb_base: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let nucs = [b'A', b'C', b'T', b'G', b'a', b'c', b't', b'g'];

    let mut writer = Vec::new();

    for i in 0..nb_reads {
        let mut record = noodles_fastx::Record::default();
        record
            .name_mut()
            .extend(format!("{}", i).as_bytes().to_vec());
        *record.description_mut() = Some(format!("random read number {}", i).as_bytes().to_vec());
        record
            .sequence_mut()
            .extend((0..nb_base).map(|_| *nucs.choose(&mut rng).unwrap()));

        record.to_writer(&mut writer).unwrap();
    }

    return writer;
}

fn generate_fastq(nb_reads: usize, nb_base: usize) -> Vec<u8> {
    let mut rng = rand::thread_rng();

    let nucs = [b'A', b'C', b'T', b'G', b'a', b'c', b't', b'g'];
    let qual: Vec<u8> = (33..126).collect();

    let mut writer = Vec::new();

    for i in 0..nb_reads {
        let mut record = noodles_fastx::Record::default();
        record
            .name_mut()
            .extend(format!("{}", i).as_bytes().to_vec());
        *record.description_mut() = Some(format!("random read number {}", i).as_bytes().to_vec());
        *record.second_description_mut() =
            Some(format!("random read number {}", i).as_bytes().to_vec());
        record
            .sequence_mut()
            .extend((0..nb_base).map(|_| *nucs.choose(&mut rng).unwrap()));
        *record.quality_mut() = Some(
            (0..nb_base)
                .map(|_| *qual.choose(&mut rng).unwrap())
                .collect(),
        );

        record.to_writer(&mut writer).unwrap();
    }

    return writer;
}

fn fasta(c: &mut Criterion) {
    let mut g = c.benchmark_group("fasta");

    let input = generate_fasta(100_000, 150);

    g.bench_function("fastx", |b| {
        b.iter(|| {
            let mut reader = noodles_fastx::Reader::new(&input[..]);
            for result in reader.records() {
                let record = result.unwrap();
                black_box(record);
            }
        })
    });

    g.bench_function("fasta", |b| {
        b.iter(|| {
            let mut reader = noodles_fasta::Reader::new(&input[..]);
            for result in reader.records() {
                let record = result.unwrap();
                black_box(record);
            }
        })
    });
}

fn fastq(c: &mut Criterion) {
    let mut g = c.benchmark_group("fastq");

    let input = generate_fastq(100_000, 150);

    g.bench_function("fastx", |b| {
        b.iter(|| {
            let mut reader = noodles_fastx::Reader::new(&input[..]);
            for result in reader.records() {
                let record = result.unwrap();
                black_box(record);
            }
        })
    });

    g.bench_function("fastq", |b| {
        b.iter(|| {
            let mut reader = noodles_fastq::Reader::new(&input[..]);
            for result in reader.records() {
                let record = result.unwrap();
                black_box(record);
            }
        })
    });
}

fn setup(c: &mut Criterion) {
    fasta(c);
    fastq(c);
}

criterion_group!(benches, setup);

criterion_main!(benches);
