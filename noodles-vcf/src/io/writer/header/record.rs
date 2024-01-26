mod key;
mod value;

use std::io::{self, Write};

use self::key::write_key;
use crate::{
    header::{
        record::{
            self,
            value::{
                map::{contig::Name, AlternativeAllele, Contig, Filter, Format, Info},
                Collection, Map,
            },
        },
        FileFormat,
    },
    record::alternate_bases::allele::Symbol,
};

fn write_record<W, K, F>(writer: &mut W, key: K, f: F) -> io::Result<()>
where
    W: Write,
    K: AsRef<str>,
    F: Fn(&mut W) -> io::Result<()>,
{
    fn write_prefix<W>(writer: &mut W) -> io::Result<()>
    where
        W: Write,
    {
        const PREFIX: &[u8] = b"##";
        writer.write_all(PREFIX)
    }

    write_prefix(writer)?;
    write_key(writer, key)?;
    write_separator(writer)?;
    f(writer)?;
    Ok(())
}

pub(super) fn write_file_format<W>(writer: &mut W, file_format: FileFormat) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::FILE_FORMAT, |w| {
        value::write_file_format(w, file_format)
    })
}

pub(super) fn write_info<W>(
    writer: &mut W,
    id: &crate::record::info::field::Key,
    info: &Map<Info>,
) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::INFO, |w| {
        value::write_map(w, id, |x| value::map::write_info(x, info))
    })
}

pub(super) fn write_filter<W>(writer: &mut W, id: &str, filter: &Map<Filter>) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::FILTER, |w| {
        value::write_map(w, id, |x| value::map::write_filter(x, filter))
    })
}

pub(super) fn write_format<W>(
    writer: &mut W,
    id: &crate::record::genotypes::keys::Key,
    format: &Map<Format>,
) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::FORMAT, |w| {
        value::write_map(w, id, |x| value::map::write_format(x, format))
    })
}

pub(super) fn write_alternative_allele<W>(
    writer: &mut W,
    id: &Symbol,
    alternative_allele: &Map<AlternativeAllele>,
) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::ALTERNATIVE_ALLELE, |w| {
        value::write_map(w, id.to_string(), |x| {
            value::map::write_alternative_allele(x, alternative_allele)
        })
    })
}

pub(super) fn write_contig<W>(writer: &mut W, id: &Name, contig: &Map<Contig>) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::CONTIG, |w| {
        value::write_map(w, id, |x| value::map::write_contig(x, contig))
    })
}

pub(super) fn write_other<W>(
    writer: &mut W,
    key: &record::key::Other,
    collection: &Collection,
) -> io::Result<()>
where
    W: Write,
{
    const META: &str = "META";

    match collection {
        Collection::Unstructured(vs) => {
            for v in vs {
                write_record(writer, key, |w| value::write_string(w, v))?;
            }
        }
        Collection::Structured(maps) => {
            for (id, map) in maps {
                write_record(writer, key, |w| {
                    value::write_map(w, id, |x| {
                        if key.as_ref() == META {
                            value::map::write_meta(x, map)
                        } else {
                            value::map::write_other(x, map)
                        }
                    })
                })?;
            }
        }
    }

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'=';
    writer.write_all(&[SEPARATOR])
}
