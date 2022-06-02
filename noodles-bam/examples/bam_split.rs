//! Splits a BAM into multiple files by read group.
//!
//! Read groups are determined by the read group records in the SAM header. Each output is named
//! `out_<index>.bam` and contains records from a single read group. Records without a read group
//! are discarded.
//!
//! This is similar to the outputs of `samtools split <src>`.

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

use std::{collections::HashMap, env, fs::File, io};

type Writers = HashMap<String, bam::Writer<bgzf::Writer<File>>>;

fn build_writers(read_groups: &sam::header::ReadGroups) -> io::Result<Writers> {
    read_groups
        .values()
        .enumerate()
        .map(|(i, rg)| {
            let dst = format!("out_{}.bam", i);
            File::create(dst).map(|f| (rg.id().into(), bam::Writer::new(f)))
        })
        .collect::<Result<_, _>>()
}

fn write_headers(writers: &mut Writers, header: &sam::Header) -> io::Result<()> {
    for read_group in header.read_groups().values() {
        let id = read_group.id();

        let writer = writers.get_mut(id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid read group: {}", id),
            )
        })?;

        let mut modified_header = header.clone();

        let read_groups = modified_header.read_groups_mut();
        read_groups.clear();
        read_groups.insert(id.into(), read_group.clone());

        writer.write_header(&modified_header)?;
        writer.write_reference_sequences(modified_header.reference_sequences())?;
    }

    Ok(())
}

fn find_read_group(data: &sam::record::Data) -> io::Result<Option<String>> {
    use sam::record::data::field::{value::Type, Tag};

    match data.get(Tag::ReadGroup) {
        Some(field) => field
            .value()
            .as_str()
            .map(|s| Some(s.into()))
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("expected {:?}, got {:?}", Type::String, field.value()),
                )
            }),
        None => Ok(None),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header: sam::Header = reader.read_header()?.parse()?;
    reader.read_reference_sequences()?;

    let mut writers = build_writers(header.read_groups())?;
    write_headers(&mut writers, &header)?;

    for result in reader.records() {
        let record = result?;

        if let Some(rg) = find_read_group(record.data())? {
            let writer = writers.get_mut(&rg).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid read group: {}", rg),
                )
            })?;

            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
