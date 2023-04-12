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
        .keys()
        .enumerate()
        .map(|(i, id)| {
            let dst = format!("out_{i}.bam");

            bam::writer::Builder::default()
                .build_from_path(dst)
                .map(|writer| (id.clone(), writer))
        })
        .collect::<Result<_, _>>()
}

fn write_headers(writers: &mut Writers, header: &sam::Header) -> io::Result<()> {
    for (id, read_group) in header.read_groups() {
        let writer = writers.get_mut(id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid read group: {id}"),
            )
        })?;

        let mut modified_header = header.clone();

        let read_groups = modified_header.read_groups_mut();
        read_groups.clear();
        read_groups.insert(id.into(), read_group.clone());

        writer.write_header(&modified_header)?;
    }

    Ok(())
}

fn find_read_group(data: &sam::record::Data) -> io::Result<Option<String>> {
    use sam::record::data::field::{Tag, Type};

    match data.get(Tag::ReadGroup) {
        Some(value) => value.as_str().map(|s| Some(s.into())).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected {:?}, got {:?}", Type::String, value),
            )
        }),
        None => Ok(None),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = bam::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let mut writers = build_writers(header.read_groups())?;
    write_headers(&mut writers, &header)?;

    for result in reader.records(&header) {
        let record = result?;

        if let Some(rg) = find_read_group(record.data())? {
            let writer = writers.get_mut(&rg).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid read group: {rg}"),
                )
            })?;

            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
