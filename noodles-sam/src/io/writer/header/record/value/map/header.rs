mod version;

use std::io::{self, Write};

use self::version::write_version_field;
use super::{write_field, write_other_fields};
use crate::header::record::value::{
    map::{header::tag, Header},
    Map,
};

pub(crate) fn write_header<W>(writer: &mut W, header: &Map<Header>) -> io::Result<()>
where
    W: Write,
{
    write_version_field(writer, header.version())?;

    if let Some(sort_order) = header.sort_order() {
        write_field(writer, tag::SORT_ORDER, sort_order)?;
    }

    if let Some(group_order) = header.group_order() {
        write_field(writer, tag::GROUP_ORDER, group_order)?;
    }

    if let Some(subsort_order) = header.subsort_order() {
        write_field(writer, tag::SUBSORT_ORDER, subsort_order)?;
    }

    write_other_fields(writer, header.other_fields())?;

    Ok(())
}
