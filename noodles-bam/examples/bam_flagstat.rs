//! Prints flag statistics of a BAM file.
//!
//! The results match the output of `samtools flagstat <src>`.

use std::{env, fmt, fs::File, io};

use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::Record};

#[derive(Debug, Default)]
struct Counts {
    read: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicate: u64,
    primary_duplicate: u64,
    mapped: u64,
    primary_mapped: u64,
    paired: u64,
    read_1: u64,
    read_2: u64,
    proper_pair: u64,
    mate_mapped: u64,
    singleton: u64,
    mate_reference_sequence_id_mismatch: u64,
    mate_reference_sequence_id_mismatch_hq: u64,
}

fn count(counts: &mut Counts, record: &Record) {
    let flags = record.flags();

    counts.read += 1;

    if !flags.is_unmapped() {
        counts.mapped += 1;
    }

    if flags.is_duplicate() {
        counts.duplicate += 1;
    }

    if flags.is_secondary() {
        counts.secondary += 1;
    } else if flags.is_supplementary() {
        counts.supplementary += 1;
    } else {
        counts.primary += 1;

        if !flags.is_unmapped() {
            counts.primary_mapped += 1;
        }

        if flags.is_duplicate() {
            counts.primary_duplicate += 1;
        }

        if flags.is_segmented() {
            counts.paired += 1;

            if flags.is_first_segment() {
                counts.read_1 += 1;
            }

            if flags.is_last_segment() {
                counts.read_2 += 1;
            }

            if !flags.is_unmapped() {
                if flags.is_properly_aligned() {
                    counts.proper_pair += 1;
                }

                if flags.is_mate_unmapped() {
                    counts.singleton += 1;
                } else {
                    counts.mate_mapped += 1;

                    if record.mate_reference_sequence_id() != record.reference_sequence_id() {
                        counts.mate_reference_sequence_id_mismatch += 1;

                        let mapq = record
                            .mapping_quality()
                            .map(u8::from)
                            .unwrap_or(sam::record::mapping_quality::MISSING);

                        if mapq >= 5 {
                            counts.mate_reference_sequence_id_mismatch_hq += 1;
                        }
                    }
                }
            }
        }
    }
}

struct PercentageFormat(u64, u64);

impl fmt::Display for PercentageFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.1 == 0 {
            f.write_str("N/A")
        } else {
            let (a, b) = (self.0 as f64, self.1 as f64);
            write!(f, "{:.2}%", a / b * 100.0)
        }
    }
}

fn print_stats(qc_pass_counts: &Counts, qc_fail_counts: &Counts) {
    println!(
        "{} + {} in total (QC-passed reads + QC-failed reads)",
        qc_pass_counts.read, qc_fail_counts.read
    );
    println!(
        "{} + {} primary",
        qc_pass_counts.primary, qc_fail_counts.primary
    );
    println!(
        "{} + {} secondary",
        qc_pass_counts.secondary, qc_fail_counts.secondary
    );
    println!(
        "{} + {} supplementary",
        qc_pass_counts.supplementary, qc_fail_counts.supplementary
    );
    println!(
        "{} + {} duplicates",
        qc_pass_counts.duplicate, qc_fail_counts.duplicate
    );
    println!(
        "{} + {} primary duplicates",
        qc_pass_counts.primary_duplicate, qc_fail_counts.primary_duplicate
    );
    println!(
        "{} + {} mapped ({} : {})",
        qc_pass_counts.mapped,
        qc_fail_counts.mapped,
        PercentageFormat(qc_pass_counts.mapped, qc_pass_counts.read),
        PercentageFormat(qc_fail_counts.mapped, qc_fail_counts.read)
    );
    println!(
        "{} + {} primary mapped ({} : {})",
        qc_pass_counts.primary_mapped,
        qc_fail_counts.primary_mapped,
        PercentageFormat(qc_pass_counts.primary_mapped, qc_pass_counts.primary),
        PercentageFormat(qc_fail_counts.primary_mapped, qc_fail_counts.primary)
    );
    println!(
        "{} + {} paired in sequencing",
        qc_pass_counts.paired, qc_fail_counts.paired
    );
    println!(
        "{} + {} read1",
        qc_pass_counts.read_1, qc_fail_counts.read_1
    );
    println!(
        "{} + {} read2",
        qc_pass_counts.read_2, qc_fail_counts.read_2
    );
    println!(
        "{} + {} properly paired ({} : {})",
        qc_pass_counts.proper_pair,
        qc_fail_counts.proper_pair,
        PercentageFormat(qc_pass_counts.proper_pair, qc_pass_counts.paired),
        PercentageFormat(qc_fail_counts.proper_pair, qc_fail_counts.paired)
    );
    println!(
        "{} + {} with itself and mate mapped",
        qc_pass_counts.mate_mapped, qc_fail_counts.mate_mapped
    );
    println!(
        "{} + {} singletons ({} : {})",
        qc_pass_counts.singleton,
        qc_fail_counts.singleton,
        PercentageFormat(qc_pass_counts.singleton, qc_pass_counts.paired),
        PercentageFormat(qc_fail_counts.singleton, qc_fail_counts.paired)
    );
    println!(
        "{} + {} with mate mapped to a different chr",
        qc_pass_counts.mate_reference_sequence_id_mismatch,
        qc_fail_counts.mate_reference_sequence_id_mismatch
    );
    println!(
        "{} + {} with mate mapped to a different chr (mapQ>=5)",
        qc_pass_counts.mate_reference_sequence_id_mismatch_hq,
        qc_fail_counts.mate_reference_sequence_id_mismatch_hq
    );
}

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    reader.read_header()?;
    reader.read_reference_sequences()?;

    let mut qc_pass_counts = Counts::default();
    let mut qc_fail_counts = Counts::default();

    for result in reader.records() {
        let record = result?;

        if record.flags().is_qc_fail() {
            count(&mut qc_fail_counts, &record);
        } else {
            count(&mut qc_pass_counts, &record);
        }
    }

    print_stats(&qc_pass_counts, &qc_fail_counts);

    Ok(())
}
