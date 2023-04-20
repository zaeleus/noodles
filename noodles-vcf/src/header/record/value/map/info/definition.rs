//! VCF header info key.

mod v4_3;
mod v4_4;

use crate::{
    header::{record::value::map::info::Type, FileFormat, Number},
    record::info::field::Key,
};

pub(crate) fn definition(
    file_format: FileFormat,
    key: &Key,
) -> Option<(Number, Type, &'static str)> {
    match key {
        Key::Standard(k) => match (file_format.major(), file_format.minor()) {
            (4, 4) => v4_4::definition(*k),
            (4, 3) => v4_3::definition(*k),
            _ => None,
        },
        Key::Other(_) => None,
    }
}
