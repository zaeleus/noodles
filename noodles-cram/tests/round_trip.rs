//! Round-trip integration tests: write CRAM records then read them back.

use std::io::Cursor;

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::io::Write as _};

mod common;

fn round_trip_cram(
    header: &sam::Header,
    records: &[sam::Record],
    reference_sequence_repository: &fasta::Repository,
    version: cram::file_definition::Version,
) -> Vec<String> {
    let cram_buf = common::write_cram(header, records, reference_sequence_repository, version);

    let cursor = Cursor::new(&cram_buf);
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(reference_sequence_repository.clone())
        .build_from_reader(cursor);

    let read_header = reader.read_header().unwrap();

    let mut sam_buf = Vec::new();
    {
        let mut sam_writer = sam::io::Writer::new(std::io::BufWriter::new(&mut sam_buf));

        for result in reader.records(&read_header) {
            let record = result.unwrap();
            sam_writer
                .write_alignment_record(&read_header, &record)
                .unwrap();
        }
    }

    String::from_utf8(sam_buf)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

fn assert_round_trip(
    header: &sam::Header,
    records: &[sam::Record],
    repository: &fasta::Repository,
    version: cram::file_definition::Version,
    label: &str,
) {
    let actual = round_trip_cram(header, records, repository, version);
    let expected = common::records_to_sam_lines(header, records);

    assert_eq!(
        actual.len(),
        expected.len(),
        "{label}: record count mismatch (got {}, expected {})",
        actual.len(),
        expected.len()
    );

    for (i, (a, e)) in actual.iter().zip(&expected).enumerate() {
        assert_eq!(
            common::normalize_sam_line(a),
            common::normalize_sam_line(e),
            "{label}: record {i} mismatch"
        );
    }
}

// ---------------------------------------------------------------------------
// Basic round-trip tests (one per version)
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_cram_v2_1() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();
    let version = cram::file_definition::Version::new(2, 1);
    assert_round_trip(&header, &records, &repository, version, "CRAM 2.1");
}

#[test]
fn test_round_trip_cram_v3_0() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();
    let version = cram::file_definition::Version::new(3, 0);
    assert_round_trip(&header, &records, &repository, version, "CRAM 3.0");
}

#[test]
fn test_round_trip_cram_v3_1() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();
    let version = cram::file_definition::Version::new(3, 1);
    assert_round_trip(&header, &records, &repository, version, "CRAM 3.1");
}

#[test]
fn test_round_trip_cram_v4_0() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();
    let version = cram::file_definition::Version::new(4, 0);
    assert_round_trip(&header, &records, &repository, version, "CRAM 4.0");
}

// ---------------------------------------------------------------------------
// Unmapped reads
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_unmapped_reads() {
    let sam_text = format!(
        "{}read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         read3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} unmapped"),
        );
    }
}

// ---------------------------------------------------------------------------
// Mixed mapped and unmapped reads
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_mixed_mapped_unmapped() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         mapped2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n\
         unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGGGG\tAAAAAAAA\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} mixed mapped/unmapped"),
        );
    }
}

// ---------------------------------------------------------------------------
// Tag variety: different tag types
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_diverse_tags() {
    let sam_text = format!(
        "{}tags1\t0\tCHROMOSOME_I\t100\t30\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXI:i:42\tXF:f:3.14\tXZ:Z:hello\tXA:A:Q\n\
         tags2\t0\tCHROMOSOME_I\t200\t30\t10M\t*\t0\t0\tGGCCAAGGTT\t??????????\tXI:i:-7\tXZ:Z:world\tXA:A:R\n\
         tags3\t0\tCHROMOSOME_I\t300\t30\t10M\t*\t0\t0\tTTTTAAAACC\t!!!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} diverse tags"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: single record
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_single_record() {
    let sam_text = format!(
        "{}single\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} single record"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: no quality scores (0xFF / *)
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_missing_quality_scores() {
    let sam_text = format!(
        "{}noqc1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\t*\n\
         noqc2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tGGCCAAGG\t*\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} missing quality scores"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: reads with soft and hard clips
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_clipped_reads() {
    let sam_text = format!(
        "{}softclip\t0\tCHROMOSOME_I\t100\t30\t3S5M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         hardclip\t0\tCHROMOSOME_I\t200\t30\t2H6M\t*\t0\t0\tGCCAAG\t??????\n\
         bothclip\t0\tCHROMOSOME_I\t300\t30\t1H2S4M2S1H\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} clipped reads"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: reads with insertions and deletions
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_indels() {
    let sam_text = format!(
        "{}insertion\t0\tCHROMOSOME_I\t100\t30\t3M2I5M\t*\t0\t0\tACGTTACGTA\tIIIIIIIIII\n\
         deletion\t0\tCHROMOSOME_I\t200\t30\t3M2D5M\t*\t0\t0\tACGACGTA\tIIIIIIII\n\
         complex\t0\tCHROMOSOME_I\t300\t30\t2M1I3M1D2M\t*\t0\t0\tACTGTAAC\tIIIIIIII\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} indels"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: reverse strand reads
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_reverse_strand() {
    let sam_text = format!(
        "{}forward\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         reverse\t16\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} reverse strand"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: empty file (header only, no records)
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_empty_file() {
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let mut cram_buf = Vec::new();
        {
            let header = sam::Header::default();
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_version(version)
                .build_from_writer(&mut cram_buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();
            writer.try_finish(&header).unwrap();
        }

        let cursor = Cursor::new(&cram_buf);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);

        let read_header = reader.read_header().unwrap();
        let record_count = reader.records(&read_header).count();

        assert_eq!(record_count, 0, "{label} empty file: expected 0 records");
    }
}

// ---------------------------------------------------------------------------
// Edge case: longer reads (150bp, typical Illumina)
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_long_reads() {
    let seq = "ACGTACGT".repeat(19); // 152bp
    let qual = "I".repeat(152);
    let sam_text = format!(
        "{}long1\t0\tCHROMOSOME_I\t100\t30\t152M\t*\t0\t0\t{seq}\t{qual}\n\
         long2\t0\tCHROMOSOME_I\t500\t30\t152M\t*\t0\t0\t{seq}\t{qual}\n",
        common::SINGLE_CHROM_HEADER,
    );

    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} long reads"),
        );
    }
}

// ---------------------------------------------------------------------------
// Edge case: multiple records on different references (multi-ref container)
// ---------------------------------------------------------------------------

#[test]
fn test_round_trip_multi_reference() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();

    // The ce#5b SAM data has reads on CHROMOSOME_I through CHROMOSOME_V,
    // so this is inherently a multi-reference test.
    for (version, label) in common::all_versions() {
        assert_round_trip(
            &header,
            &records,
            &repository,
            version,
            &format!("{label} multi-reference"),
        );
    }
}

// ---------------------------------------------------------------------------
// Cross-compatibility: write with noodles, read with samtools
// ---------------------------------------------------------------------------

#[test]
#[ignore] // Requires samtools to be installed locally
fn test_samtools_can_read_noodles_cram() {
    use std::io::Write;
    use std::process::Command;

    // Check if samtools is available
    let samtools = Command::new("samtools").arg("--version").output();
    if samtools.is_err() || !samtools.unwrap().status.success() {
        eprintln!("samtools not found, skipping cross-compatibility test");
        return;
    }

    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();

    // Write a temp FASTA file for samtools
    let tmp_dir = std::env::temp_dir();
    let pid = std::process::id();
    let fasta_path = tmp_dir.join(format!("noodles_test_ref_{pid}.fa"));
    {
        let mut fa = std::io::BufWriter::new(std::fs::File::create(&fasta_path).unwrap());
        for (name, len) in [
            ("CHROMOSOME_I", 1000),
            ("CHROMOSOME_II", 200),
            ("CHROMOSOME_III", 200),
            ("CHROMOSOME_IV", 200),
            ("CHROMOSOME_V", 200),
        ] {
            let seq = repository.get(name.as_bytes()).unwrap().unwrap();
            writeln!(fa, ">{name}").unwrap();
            let bytes: &[u8] = (*seq).as_ref();
            for chunk in bytes[..len].chunks(80) {
                fa.write_all(chunk).unwrap();
                fa.write_all(b"\n").unwrap();
            }
        }
    }

    // Index the FASTA for samtools
    let idx_result = Command::new("samtools")
        .args(["faidx", &fasta_path.to_string_lossy()])
        .output()
        .unwrap();
    assert!(
        idx_result.status.success(),
        "samtools faidx failed: {}",
        String::from_utf8_lossy(&idx_result.stderr)
    );

    let versions = common::all_versions();

    for (version, label) in versions {
        // Write CRAM with noodles
        let cram_buf = common::write_cram(&header, &records, &repository, version);

        // Write CRAM to a temp file
        let cram_path = tmp_dir.join(format!(
            "noodles_test_{pid}_{}.cram",
            label.replace(' ', "_")
        ));
        std::fs::write(&cram_path, &cram_buf).unwrap();

        // Read with samtools view
        let output = Command::new("samtools")
            .args(["view", "-h", "--no-PG"])
            .arg("--reference")
            .arg(&fasta_path)
            .arg(&cram_path)
            .output()
            .unwrap();

        std::fs::remove_file(&cram_path).ok();

        assert!(
            output.status.success(),
            "{label}: samtools failed to read noodles CRAM: {}",
            String::from_utf8_lossy(&output.stderr)
        );

        let samtools_output = String::from_utf8(output.stdout).unwrap();
        let samtools_records: Vec<&str> = samtools_output
            .lines()
            .filter(|l| !l.starts_with('@'))
            .collect();

        assert_eq!(
            samtools_records.len(),
            records.len(),
            "{label}: samtools record count mismatch"
        );

        let expected = common::records_to_sam_lines(&header, &records);

        for (i, (actual_line, expected_line)) in samtools_records.iter().zip(&expected).enumerate()
        {
            assert_eq!(
                common::normalize_sam_line(actual_line),
                common::normalize_sam_line(expected_line),
                "{label}: samtools record {i} mismatch"
            );
        }
    }

    // Cleanup
    std::fs::remove_file(&fasta_path).ok();
    std::fs::remove_file(fasta_path.with_extension("fa.fai")).ok();
}
