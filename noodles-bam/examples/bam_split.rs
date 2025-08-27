//! Splits a BAM into multiple files by read group.
//!
//! Read groups are determined by the read group records in the SAM header. Each output is named
//! `out_<index>.bam` and contains records from a single read group. Records without a read group
//! are discarded.
//!
//! This is similar to the outputs of `samtools split <src>`.

use bstr::{BStr, ByteSlice};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

use std::{collections::HashMap, env, fs::File, io};

type Writers<'h> = HashMap<&'h BStr, bam::io::Writer<bgzf::io::Writer<File>>>;

fn build_writers(read_groups: &sam::header::ReadGroups) -> io::Result<Writers<'_>> {
    read_groups
        .keys()
        .enumerate()
        .map(|(i, id)| {
            let dst = format!("out_{i}.bam");

            File::create(dst)
                .map(bam::io::Writer::new)
                .map(|writer| (id.as_ref(), writer))
        })
        .collect::<Result<_, _>>()
}

fn write_headers(writers: &mut Writers<'_>, header: &sam::Header) -> io::Result<()> {
    for (id, read_group) in header.read_groups() {
        let writer = writers.get_mut(id.as_bstr()).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid read group: {id}"),
            )
        })?;

        let mut modified_header = header.clone();

        let read_groups = modified_header.read_groups_mut();
        read_groups.clear();
        read_groups.insert(id.clone(), read_group.clone());

        writer.write_header(&modified_header)?;
    }

    Ok(())
}

fn get_read_group(data: &dyn sam::alignment::record::Data) -> Option<io::Result<&BStr>> {
    use sam::alignment::record::data::field::{Tag, Type, Value};

    data.get(&Tag::READ_GROUP).map(|result| {
        result.and_then(|value| match value {
            Value::String(s) => Ok(s),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected {:?}, got {:?}", Type::String, value.ty()),
            )),
        })
    })
}

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let mut writers = build_writers(header.read_groups())?;
    write_headers(&mut writers, &header)?;

    for result in reader.records() {
        let record = result?;

        if let Some(read_group) = get_read_group(&record.data()).transpose()? {
            let writer = writers.get_mut(read_group).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid read group: {read_group:?}"),
                )
            })?;

            writer.write_record(&header, &record)?;
        }
    }

    Ok(())
}
