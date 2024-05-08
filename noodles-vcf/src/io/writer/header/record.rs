mod key;
mod value;

use std::io::{self, Write};

use self::key::write_key;
use super::write_newline;
use crate::header::{
    record::{
        self,
        value::{
            map::{AlternativeAllele, Contig, Filter, Format, Info},
            Collection, Map,
        },
    },
    FileFormat,
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
    write_newline(writer)?;

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

pub(super) fn write_info<W>(writer: &mut W, id: &str, info: &Map<Info>) -> io::Result<()>
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

pub(super) fn write_format<W>(writer: &mut W, id: &str, format: &Map<Format>) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::FORMAT, |w| {
        value::write_map(w, id, |x| value::map::write_format(x, format))
    })
}

pub(super) fn write_alternative_allele<W>(
    writer: &mut W,
    id: &str,
    alternative_allele: &Map<AlternativeAllele>,
) -> io::Result<()>
where
    W: Write,
{
    write_record(writer, &record::key::ALTERNATIVE_ALLELE, |w| {
        value::write_map(w, id, |x| {
            value::map::write_alternative_allele(x, alternative_allele)
        })
    })
}

pub(super) fn write_contig<W>(writer: &mut W, id: &str, contig: &Map<Contig>) -> io::Result<()>
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
                    value::write_other_map(w, map.id_tag(), id, |x| {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_other() -> Result<(), Box<dyn std::error::Error>> {
        let mut buf = Vec::new();

        let key = "comment".parse()?;

        buf.clear();
        let collection =
            Collection::Unstructured(vec![String::from("noodles"), String::from("vcf")]);
        write_other(&mut buf, &key, &collection)?;
        assert_eq!(buf, b"##comment=noodles\n##comment=vcf\n");

        buf.clear();
        let collection = Collection::Structured(
            [
                (String::from("noodles"), Map::default()),
                (String::from("vcf"), Map::default()),
            ]
            .into_iter()
            .collect(),
        );
        write_other(&mut buf, &key, &collection)?;
        assert_eq!(buf, b"##comment=<ID=noodles>\n##comment=<ID=vcf>\n");

        Ok(())
    }
}
