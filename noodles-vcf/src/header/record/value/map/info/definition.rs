//! VCF header info key.

mod v4_3;
mod v4_4;
mod v4_5;

use crate::header::{
    record::value::map::info::{Number, Type},
    FileFormat,
};

pub(crate) fn definition(
    file_format: FileFormat,
    key: &str,
) -> Option<(Number, Type, &'static str)> {
    match (file_format.major(), file_format.minor()) {
        (4, 5) => v4_5::definition(key),
        (4, 4) => v4_4::definition(key),
        (4, 3) => v4_3::definition(key),
        _ => None,
    }
}
