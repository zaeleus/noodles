//! Integration tests for Reader/IndexedReader methods and writer builder options.
//!
//! Covers:
//!   - `Reader::for_each_record()`
//!   - `Reader::query_unmapped()`
//!   - `IndexedReader` delegation of the above
//!   - Writer builder options (strip_md_nm, embed_reference_sequences, etc.)
//!   - `Reader::query()` region-based queries

use std::{
    io::{self, BufWriter, Cursor},
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};

use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::{self as sam, alignment::io::Write as _};

mod common;

// ---------------------------------------------------------------------------
// Test-specific helpers (not shared)
// ---------------------------------------------------------------------------

static COUNTER: AtomicUsize = AtomicUsize::new(0);

/// Generate a unique temp file path.
fn unique_tmp_path(prefix: &str) -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::SeqCst);
    let pid = std::process::id();
    std::env::temp_dir().join(format!("{prefix}_{pid}_{id}.cram"))
}

/// Build a CRAI index from in-memory CRAM data by writing to a temp file,
/// running `cram::fs::index()`, then cleaning up.
fn build_index(cram_data: &[u8]) -> cram::crai::Index {
    let path = unique_tmp_path("noodles_test_index");
    std::fs::write(&path, cram_data).unwrap();
    let index = cram::fs::index(&path).unwrap();
    std::fs::remove_file(&path).ok();
    index
}

// ---------------------------------------------------------------------------
// Test Group 1: for_each_record()
// ---------------------------------------------------------------------------

#[test]
fn test_for_each_record_matches_records_iterator() {
    let (header, records) = common::parse_sam(common::CE5B_SAM);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        // Collect via records() iterator
        let iter_lines: Vec<String> = {
            let cursor = Cursor::new(&cram_data);
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .build_from_reader(cursor);
            let read_header = reader.read_header().unwrap();
            let mut buf = Vec::new();
            {
                let mut w = sam::io::Writer::new(BufWriter::new(&mut buf));
                for result in reader.records(&read_header) {
                    let record = result.unwrap();
                    w.write_alignment_record(&read_header, &record).unwrap();
                }
            }
            String::from_utf8(buf)
                .unwrap()
                .lines()
                .map(String::from)
                .collect()
        };

        // Collect via for_each_record()
        let foreach_lines: Vec<String> = {
            let cursor = Cursor::new(&cram_data);
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .build_from_reader(cursor);
            let read_header = reader.read_header().unwrap();
            let mut buf = Vec::new();
            {
                let mut w = sam::io::Writer::new(BufWriter::new(&mut buf));
                reader
                    .for_each_record(&read_header, |record| {
                        w.write_alignment_record(&read_header, record)
                    })
                    .unwrap();
            }
            String::from_utf8(buf)
                .unwrap()
                .lines()
                .map(String::from)
                .collect()
        };

        assert_eq!(
            iter_lines.len(),
            foreach_lines.len(),
            "{label}: record count mismatch between records() and for_each_record()"
        );

        for (i, (a, b)) in iter_lines.iter().zip(&foreach_lines).enumerate() {
            assert_eq!(
                common::normalize_sam_line(a),
                common::normalize_sam_line(b),
                "{label}: record {i} mismatch"
            );
        }
    }
}

#[test]
fn test_for_each_record_empty_file() {
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let header = sam::Header::default();
        let cram_data = common::write_cram(&header, &[], &repository, version);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut count = 0usize;
        reader
            .for_each_record(&read_header, |_record| {
                count += 1;
                Ok(())
            })
            .unwrap();

        assert_eq!(count, 0, "{label}: expected 0 records in empty file");
    }
}

#[test]
fn test_for_each_record_unmapped_only() {
    let sam_text = format!(
        "{}read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         read3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut count = 0usize;
        let mut all_unmapped = true;
        reader
            .for_each_record(&read_header, |record| {
                count += 1;
                let flags = record.flags()?;
                if !flags.is_unmapped() {
                    all_unmapped = false;
                }
                Ok(())
            })
            .unwrap();

        assert_eq!(count, 3, "{label}: expected 3 unmapped records");
        assert!(all_unmapped, "{label}: all records should be unmapped");
    }
}

#[test]
fn test_for_each_record_mixed_mapped_unmapped() {
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
        let cram_data = common::write_cram(&header, &records, &repository, version);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut mapped_count = 0usize;
        let mut unmapped_count = 0usize;
        reader
            .for_each_record(&read_header, |record| {
                let flags = record.flags()?;
                if flags.is_unmapped() {
                    unmapped_count += 1;
                } else {
                    mapped_count += 1;
                }
                Ok(())
            })
            .unwrap();

        assert_eq!(mapped_count, 2, "{label}: expected 2 mapped records");
        assert_eq!(unmapped_count, 2, "{label}: expected 2 unmapped records");
    }
}

#[test]
fn test_for_each_record_closure_error_propagation() {
    let sam_text = format!(
        "{}rec1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         rec2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tGGCCAAGG\t????????\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut call_count = 0usize;
        let result = reader.for_each_record(&read_header, |_record| {
            call_count += 1;
            Err(io::Error::other("test error"))
        });

        assert!(result.is_err(), "{label}: expected error from closure");
        assert_eq!(
            result.unwrap_err().to_string(),
            "test error",
            "{label}: error message mismatch"
        );
        assert_eq!(
            call_count, 1,
            "{label}: closure should be called once before error"
        );
    }
}

#[test]
fn test_for_each_record_record_data_accessible() {
    let sam_text = format!(
        "{}testread\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut saw_record = false;
        reader
            .for_each_record(&read_header, |record| {
                saw_record = true;

                // Check name
                let name = record.name().expect("record should have a name");
                assert_eq!(name, b"testread".as_slice(), "{label}: name mismatch");

                // Check flags (not unmapped)
                let flags = record.flags()?;
                assert!(!flags.is_unmapped(), "{label}: should be mapped");

                // Check alignment_start
                let start = record
                    .alignment_start()
                    .expect("should have alignment start")?;
                assert_eq!(usize::from(start), 100, "{label}: alignment start mismatch");

                Ok(())
            })
            .unwrap();

        assert!(saw_record, "{label}: closure should have been called");
    }
}

// ---------------------------------------------------------------------------
// Test Group 2: query_unmapped()
// ---------------------------------------------------------------------------

#[test]
fn test_query_unmapped_mixed() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         mapped2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGGGG\tAAAAAAAA\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(unmapped.len(), 2, "{label}: expected 2 unmapped records");
        for rec in &unmapped {
            assert!(
                rec.flags().is_unmapped(),
                "{label}: record should be unmapped"
            );
        }
    }
}

#[test]
fn test_query_unmapped_all_unmapped() {
    let sam_text = format!(
        "{}read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         read3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(unmapped.len(), 3, "{label}: expected 3 unmapped records");
    }
}

#[test]
fn test_query_unmapped_no_unmapped() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         mapped2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(unmapped.len(), 0, "{label}: expected 0 unmapped records");
    }
}

#[test]
fn test_query_unmapped_empty_file() {
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let header = sam::Header::default();
        let cram_data = common::write_cram(&header, &[], &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(unmapped.len(), 0, "{label}: expected 0 unmapped records");
    }
}

#[test]
fn test_query_unmapped_verifies_record_content() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         uname_a\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         uname_b\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGGGG\tAAAAAAAA\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(unmapped.len(), 2, "{label}: expected 2 unmapped records");

        let mut names: Vec<String> = unmapped
            .iter()
            .map(|r| {
                r.name()
                    .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
                    .unwrap_or_default()
            })
            .collect();
        names.sort();

        assert_eq!(
            names,
            vec!["uname_a", "uname_b"],
            "{label}: unmapped record names mismatch"
        );
    }
}

#[test]
fn test_query_unmapped_with_empty_index() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        // Use an empty index — no unmapped entries, so seeks to EOF
        let empty_index = cram::crai::Index::default();

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &empty_index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(
            unmapped.len(),
            0,
            "{label}: expected 0 records with empty index"
        );
    }
}

// ---------------------------------------------------------------------------
// Test Group 3: IndexedReader delegation
// ---------------------------------------------------------------------------

#[test]
fn test_indexed_reader_for_each_record() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         mapped2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .set_index(index)
            .build_from_reader(cursor)
            .unwrap();
        let read_header = reader.read_header().unwrap();

        let mut count = 0usize;
        reader
            .for_each_record(&read_header, |_record| {
                count += 1;
                Ok(())
            })
            .unwrap();

        assert_eq!(
            count, 3,
            "{label}: expected 3 records via IndexedReader::for_each_record"
        );
    }
}

#[test]
fn test_indexed_reader_query_unmapped() {
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n\
         unmapped3\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGGGG\tAAAAAAAA\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .set_index(index)
            .build_from_reader(cursor)
            .unwrap();
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(
            unmapped.len(),
            3,
            "{label}: expected 3 unmapped records via IndexedReader"
        );
        for rec in &unmapped {
            assert!(
                rec.flags().is_unmapped(),
                "{label}: record should be unmapped"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test Group 4: Writer builder options
// ---------------------------------------------------------------------------

#[test]
fn test_strip_md_nm() {
    use sam::alignment::record::data::field::Tag;

    // Write records that contain MD/NM tags
    let sam_text = format!(
        "{}rec1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tMD:Z:8\tNM:i:0\n\
         rec2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\tMD:Z:8\tNM:i:0\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_version(version)
                .strip_md_nm(true)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        // Read back and verify MD/NM tags are absent
        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        for (i, result) in reader.records(&read_header).enumerate() {
            let record = result.unwrap();
            let data = record.data();

            assert!(
                data.get(&Tag::MISMATCHED_POSITIONS).is_none(),
                "{label}: record {i} should not have MD tag"
            );
            assert!(
                data.get(&Tag::EDIT_DISTANCE).is_none(),
                "{label}: record {i} should not have NM tag"
            );
        }
    }
}

#[test]
fn test_embed_reference_sequences() {
    let sam_text = format!(
        "{}rec1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         rec2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();
    let expected = common::records_to_sam_lines(&header, &records);

    for (version, label) in common::all_versions() {
        // Write with embed_reference_sequences(true) using the real repository
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_version(version)
                .embed_reference_sequences(true)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        // Read back with an empty repository — should still decode because ref is embedded
        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(fasta::Repository::default())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut actual = Vec::new();
        {
            let mut w = sam::io::Writer::new(BufWriter::new(&mut actual));
            for result in reader.records(&read_header) {
                let record = result.unwrap();
                w.write_alignment_record(&read_header, &record).unwrap();
            }
        }

        let actual_lines: Vec<String> = String::from_utf8(actual)
            .unwrap()
            .lines()
            .map(String::from)
            .collect();

        assert_eq!(
            actual_lines.len(),
            expected.len(),
            "{label}: record count mismatch for embed_reference_sequences"
        );

        for (i, (a, e)) in actual_lines.iter().zip(&expected).enumerate() {
            assert_eq!(
                common::normalize_sam_line(a),
                common::normalize_sam_line(e),
                "{label}: record {i} mismatch for embed_reference_sequences"
            );
        }
    }
}

#[test]
fn test_reference_free_mode() {
    let sam_text = format!(
        "{}read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         read3\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let expected = common::records_to_sam_lines(&header, &records);

    for (version, label) in common::all_versions() {
        // Write with no reference repository and reference_required=false
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_required(false)
                .set_version(version)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        // Read back — also without a reference
        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default().build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut actual = Vec::new();
        {
            let mut w = sam::io::Writer::new(BufWriter::new(&mut actual));
            for result in reader.records(&read_header) {
                let record = result.unwrap();
                w.write_alignment_record(&read_header, &record).unwrap();
            }
        }

        let actual_lines: Vec<String> = String::from_utf8(actual)
            .unwrap()
            .lines()
            .map(String::from)
            .collect();

        assert_eq!(
            actual_lines.len(),
            expected.len(),
            "{label}: record count mismatch for reference_free mode"
        );

        for (i, (a, e)) in actual_lines.iter().zip(&expected).enumerate() {
            assert_eq!(
                common::normalize_sam_line(a),
                common::normalize_sam_line(e),
                "{label}: record {i} mismatch for reference_free mode"
            );
        }
    }
}

#[test]
fn test_records_per_slice() {
    // Write 20 records with records_per_slice=5
    let mut sam_lines = common::SINGLE_CHROM_HEADER.to_string();
    for i in 0..20 {
        let pos = 100 + i * 10;
        sam_lines.push_str(&format!(
            "rec{i}\t0\tCHROMOSOME_I\t{pos}\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
        ));
    }
    let (header, records) = common::parse_sam(&sam_lines);
    let repository = common::make_reference_repository();
    let expected = common::records_to_sam_lines(&header, &records);

    for (version, label) in common::all_versions() {
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_version(version)
                .set_records_per_slice(5)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut actual = Vec::new();
        {
            let mut w = sam::io::Writer::new(BufWriter::new(&mut actual));
            for result in reader.records(&read_header) {
                let record = result.unwrap();
                w.write_alignment_record(&read_header, &record).unwrap();
            }
        }

        let actual_lines: Vec<String> = String::from_utf8(actual)
            .unwrap()
            .lines()
            .map(String::from)
            .collect();

        assert_eq!(
            actual_lines.len(),
            expected.len(),
            "{label}: record count mismatch for records_per_slice=5"
        );

        for (i, (a, e)) in actual_lines.iter().zip(&expected).enumerate() {
            assert_eq!(
                common::normalize_sam_line(a),
                common::normalize_sam_line(e),
                "{label}: record {i} mismatch for records_per_slice=5"
            );
        }
    }
}

#[test]
fn test_slices_per_container() {
    // Write 15 records with records_per_slice=5, slices_per_container=2
    let mut sam_lines = common::SINGLE_CHROM_HEADER.to_string();
    for i in 0..15 {
        let pos = 100 + i * 10;
        sam_lines.push_str(&format!(
            "rec{i}\t0\tCHROMOSOME_I\t{pos}\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"
        ));
    }
    let (header, records) = common::parse_sam(&sam_lines);
    let repository = common::make_reference_repository();
    let expected = common::records_to_sam_lines(&header, &records);

    for (version, label) in common::all_versions() {
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_version(version)
                .set_records_per_slice(5)
                .set_slices_per_container(2)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let mut actual = Vec::new();
        {
            let mut w = sam::io::Writer::new(BufWriter::new(&mut actual));
            for result in reader.records(&read_header) {
                let record = result.unwrap();
                w.write_alignment_record(&read_header, &record).unwrap();
            }
        }

        let actual_lines: Vec<String> = String::from_utf8(actual)
            .unwrap()
            .lines()
            .map(String::from)
            .collect();

        assert_eq!(
            actual_lines.len(),
            expected.len(),
            "{label}: record count mismatch for slices_per_container=2"
        );

        for (i, (a, e)) in actual_lines.iter().zip(&expected).enumerate() {
            assert_eq!(
                common::normalize_sam_line(a),
                common::normalize_sam_line(e),
                "{label}: record {i} mismatch for slices_per_container=2"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test Group 5: query() region queries
// ---------------------------------------------------------------------------

#[test]
fn test_query_region() {
    let sam_text = format!(
        "{}rec_100\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         rec_200\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tGGCCAAGG\t????????\n\
         rec_300\t0\tCHROMOSOME_I\t300\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        // Query CHROMOSOME_I:150-250 — should overlap rec_200 (200..207) and
        // possibly rec_100 (100..107, ends before 150) depending on the implementation.
        // At minimum rec_200 must be returned.
        let region: noodles_core::Region = "CHROMOSOME_I:150-250".parse().unwrap();
        let query_results: Vec<_> = reader
            .query(&read_header, &index, &region)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert!(
            !query_results.is_empty(),
            "{label}: query should return at least one record"
        );

        // All returned records should overlap the query region
        let mut found_rec_200 = false;
        for rec in &query_results {
            if let Some(name) = rec.name()
                && name == b"rec_200".as_slice()
            {
                found_rec_200 = true;
            }
        }

        assert!(
            found_rec_200,
            "{label}: query result should include rec_200"
        );
    }
}

#[test]
fn test_query_region_no_results() {
    let sam_text = format!(
        "{}rec_100\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         rec_200\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tGGCCAAGG\t????????\n\
         rec_300\t0\tCHROMOSOME_I\t300\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);
        let index = build_index(&cram_data);

        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        // Query a region with no records
        let region: noodles_core::Region = "CHROMOSOME_I:500-600".parse().unwrap();
        let query_results: Vec<_> = reader
            .query(&read_header, &index, &region)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(
            query_results.len(),
            0,
            "{label}: query should return 0 records for region with no data"
        );
    }
}

// ---------------------------------------------------------------------------
// Test Group 6: Bug fix verification
// ---------------------------------------------------------------------------

#[test]
fn test_reference_free_mode_with_sq_lines() {
    // Verifies Issue 3 fix: write_file_header no longer panics when
    // reference_required=false and the header has @SQ lines.
    let sam_text = format!(
        "{}read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         read2\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);

    for (version, label) in common::all_versions() {
        let mut buf = Vec::new();
        {
            let mut writer = cram::io::writer::Builder::default()
                .set_reference_required(false)
                .set_version(version)
                .build_from_writer(&mut buf);

            writer.write_file_definition().unwrap();
            writer.write_file_header(&header).unwrap();

            for record in &records {
                writer.write_alignment_record(&header, record).unwrap();
            }

            writer.try_finish(&header).unwrap();
        }

        // Read back and verify records round-trip
        let cursor = Cursor::new(&buf);
        let mut reader = cram::io::reader::Builder::default().build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let read_records: Vec<_> = reader
            .records(&read_header)
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(
            read_records.len(),
            2,
            "{label}: expected 2 records in reference-free mode with @SQ lines"
        );
        for rec in &read_records {
            assert!(
                rec.flags().is_unmapped(),
                "{label}: record should be unmapped"
            );
        }
    }
}

#[test]
fn test_index_mixed_data_directly() {
    // Verifies Issue 1 fix: cram::fs::index() no longer panics on CRAM files
    // with mixed mapped+unmapped records in the same container (multi-ref slices).
    let sam_text = format!(
        "{}mapped1\t0\tCHROMOSOME_I\t100\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n\
         mapped2\t0\tCHROMOSOME_I\t200\t30\t8M\t*\t0\t0\tTTTTAAAA\t!!!!!!!!\n\
         unmapped1\t4\t*\t0\t0\t*\t*\t0\t0\tGGCCAAGG\t????????\n\
         unmapped2\t4\t*\t0\t0\t*\t*\t0\t0\tCCCCGGGG\tAAAAAAAA\n",
        common::SINGLE_CHROM_HEADER,
    );
    let (header, records) = common::parse_sam(&sam_text);
    let repository = common::make_reference_repository();

    for (version, label) in common::all_versions() {
        let cram_data = common::write_cram(&header, &records, &repository, version);

        // build_index uses cram::fs::index() which reads without a reference
        // repository. Before Issue 1 fix, this would panic on multi-ref slices.
        let index = build_index(&cram_data);

        assert!(
            !index.is_empty(),
            "{label}: index should have entries for mixed data"
        );

        // Verify the index is usable for queries
        let cursor = Cursor::new(&cram_data);
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(cursor);
        let read_header = reader.read_header().unwrap();

        let unmapped: Vec<_> = reader
            .query_unmapped(&read_header, &index)
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

        assert_eq!(
            unmapped.len(),
            2,
            "{label}: expected 2 unmapped records from mixed-data index"
        );
    }
}
