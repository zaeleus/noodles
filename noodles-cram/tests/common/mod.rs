use std::io::BufWriter;

use noodles_cram as cram;
use noodles_fasta::{
    self as fasta,
    record::{Definition, Sequence},
};
use noodles_sam::{self as sam, alignment::io::Write as _};

/// Builds an in-memory reference repository with synthetic sequences matching
/// the repeating CCTAAG pattern from the original C. elegans test data.
///
/// - CHROMOSOME_I: 1000 bp (starts with G, then repeating CCTAAG)
/// - CHROMOSOME_II–V: 200 bp each (repeating CCTAAG)
pub fn make_reference_repository() -> fasta::Repository {
    fn make_sequence(prefix: &[u8], repeat: &[u8], total_len: usize) -> Vec<u8> {
        let mut seq = Vec::with_capacity(total_len);
        seq.extend_from_slice(prefix);
        while seq.len() < total_len {
            let remaining = total_len - seq.len();
            let chunk = &repeat[..remaining.min(repeat.len())];
            seq.extend_from_slice(chunk);
        }
        seq
    }

    let records = vec![
        fasta::Record::new(
            Definition::new("CHROMOSOME_I", None),
            Sequence::from(make_sequence(b"G", b"CCTAAG", 1000)),
        ),
        fasta::Record::new(
            Definition::new("CHROMOSOME_II", None),
            Sequence::from(make_sequence(b"", b"CCTAAG", 200)),
        ),
        fasta::Record::new(
            Definition::new("CHROMOSOME_III", None),
            Sequence::from(make_sequence(b"", b"CCTAAG", 200)),
        ),
        fasta::Record::new(
            Definition::new("CHROMOSOME_IV", None),
            Sequence::from(make_sequence(b"", b"CCTAAG", 200)),
        ),
        fasta::Record::new(
            Definition::new("CHROMOSOME_V", None),
            Sequence::from(make_sequence(b"", b"CCTAAG", 200)),
        ),
    ];

    fasta::Repository::new(records)
}

/// SAM header for tests using only CHROMOSOME_I.
pub const SINGLE_CHROM_HEADER: &str = "@SQ\tSN:CHROMOSOME_I\tLN:1000\n";

/// The ce#5b.sam records as inline SAM text (header + 7 records).
/// Originally from htslib test data. Records span 5 chromosomes with
/// complex CIGARs (27M1D73M, 7S20M1D23M10I30M10S).
pub const CE5B_SAM: &str = "\
@SQ\tSN:CHROMOSOME_I\tLN:1000
@SQ\tSN:CHROMOSOME_II\tLN:200
@SQ\tSN:CHROMOSOME_III\tLN:200
@SQ\tSN:CHROMOSOME_IV\tLN:200
@SQ\tSN:CHROMOSOME_V\tLN:200
I\t16\tCHROMOSOME_I\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tXS:i:-18\tAS:i:-18\tYT:Z:UU
II.14978392\t16\tCHROMOSOME_II\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tXS:i:-18\tAS:i:-18\tYT:Z:UU
III\t16\tCHROMOSOME_III\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tXS:i:-18\tAS:i:-18\tYT:Z:UU
IV\t16\tCHROMOSOME_IV\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tXS:i:-18\tAS:i:-18\tYT:Z:UU
V\t16\tCHROMOSOME_V\t2\t1\t27M1D73M\t*\t0\t0\tCCTAGCCCTAACCCTAACCCTAACCCTAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t#############################@B?8B?BA@@DDBCDDCBC@CDCDCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\tXG:i:1\tXM:i:5\tXN:i:0\tXO:i:1\tXS:i:-18\tAS:i:-18\tYT:Z:UU
VI\t0\tCHROMOSOME_V\t10\t1\t7S20M1D23M10I30M10S\t*\t0\t0\tAGCCTAAGCCTAAGCCTAAGCCTAAGCTAAGCCTAAGCCTAAGCCTAAGCTTTTTTTTTTCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAAGCCTAA\t*
VI\t256\tCHROMOSOME_V\t10\t1\t7S20M1D23M10I30M10S\t*\t0\t0\t*\t*
";

pub fn parse_sam(sam_text: &str) -> (sam::Header, Vec<sam::Record>) {
    let mut reader = sam::io::Reader::new(std::io::BufReader::new(sam_text.as_bytes()));
    let header = reader.read_header().unwrap();
    // Filter out records with empty sequences (e.g. secondary alignments with SEQ='*')
    // since CRAM round-trip cannot preserve them without sequence data.
    let records: Vec<_> = reader
        .records()
        .map(|r| r.unwrap())
        .filter(|r| !r.sequence().is_empty())
        .collect();
    (header, records)
}

/// Returns all CRAM versions supported for writing.
///
/// CRAM 2.0 is omitted because noodles always writes the minor version bump (2.1)
/// which is the de facto standard; 2.0 and 2.1 share the same on-disk format.
pub fn all_versions() -> Vec<(cram::file_definition::Version, &'static str)> {
    vec![
        (cram::file_definition::Version::new(2, 1), "CRAM 2.1"),
        (cram::file_definition::Version::new(3, 0), "CRAM 3.0"),
        (cram::file_definition::Version::new(3, 1), "CRAM 3.1"),
        (cram::file_definition::Version::new(4, 0), "CRAM 4.0"),
    ]
}

pub fn normalize_sam_line(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    let flags: u16 = fields.get(1).and_then(|f| f.parse().ok()).unwrap_or(0);
    let is_unmapped = flags & 0x4 != 0;
    let mut result = Vec::new();

    for (i, field) in fields[..11.min(fields.len())].iter().enumerate() {
        if i == 4 && is_unmapped {
            result.push("255");
            continue;
        }
        result.push(field);
    }

    for field in fields.iter().skip(11) {
        if !field.starts_with("MD:") && !field.starts_with("NM:") {
            result.push(field);
        }
    }

    result.join("\t")
}

pub fn records_to_sam_lines(header: &sam::Header, records: &[sam::Record]) -> Vec<String> {
    let mut buf = Vec::new();
    {
        let mut writer = sam::io::Writer::new(BufWriter::new(&mut buf));
        for record in records {
            writer.write_alignment_record(header, record).unwrap();
        }
    }
    String::from_utf8(buf)
        .unwrap()
        .lines()
        .map(String::from)
        .collect()
}

/// Write CRAM data to an in-memory buffer.
pub fn write_cram(
    header: &sam::Header,
    records: &[sam::Record],
    repo: &fasta::Repository,
    version: cram::file_definition::Version,
) -> Vec<u8> {
    let mut buf = Vec::new();
    {
        let mut writer = cram::io::writer::Builder::default()
            .set_reference_sequence_repository(repo.clone())
            .set_version(version)
            .build_from_writer(&mut buf);

        writer.write_file_definition().unwrap();
        writer.write_file_header(header).unwrap();

        for record in records {
            writer.write_alignment_record(header, record).unwrap();
        }

        writer.try_finish(header).unwrap();
    }
    buf
}
