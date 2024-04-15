//! Variant record samples.

pub mod keys;
mod sample;
pub mod series;

use std::io;

pub use self::{sample::Sample, series::Series};
use crate::Header;

#[allow(clippy::tabs_in_doc_comments)]
/// Variant record samples.
///
/// Variant record samples are described similarly to a data frame: rows and columns are samples
/// and series, respectively.
///
/// Take, for example, the following samples in a VCF record.
///
/// ```text
/// FORMAT	sample1	sample2	sample3
/// GT:GQ:DP	0|0:21:1	1|0:34:2	1/1:55:5
/// ```
///
/// This trait provides an interface over this data as a table.
///
/// ```text
///           GT     GQ    DP
///         +-----+-----+-----+
/// sample1 | 0|0 |  21 |   1 |
/// sample2 | 1|0 |  34 |   2 |
/// sample3 | 1/1 |  55 |   5 |
///         +-----+-----+-----+
/// ```
///
/// A sample contains the values of a row (e.g., sample1 = ["0|0", 21, 1]); and a column, the
/// values of a series (e.g., GT = ["0|0", "1|0", "1/1"]).
pub trait Samples {
    /// Returns whether there are any samples.
    fn is_empty(&self) -> bool;

    /// Returns the number of samples.
    fn len(&self) -> usize;

    /// Returns the column names.
    fn column_names<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a>;

    /// Returns the series with the given column name.
    fn select<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn Series + 'a>>>;

    /// Returns an iterator over series.
    fn series(&self) -> Box<dyn Iterator<Item = io::Result<Box<dyn Series + '_>>> + '_>;

    /// Returns an iterator over samples.
    fn iter(&self) -> Box<dyn Iterator<Item = Box<dyn Sample + '_>> + '_>;
}

impl Samples for Box<dyn Samples + '_> {
    fn is_empty(&self) -> bool {
        (**self).is_empty()
    }

    fn len(&self) -> usize {
        (**self).len()
    }

    fn column_names<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
    ) -> Box<dyn Iterator<Item = io::Result<&'a str>> + 'a> {
        (**self).column_names(header)
    }

    fn select<'a, 'h: 'a>(
        &'a self,
        header: &'h Header,
        column_name: &str,
    ) -> Option<io::Result<Box<dyn Series + 'a>>> {
        (**self).select(header, column_name)
    }

    fn series(&self) -> Box<dyn Iterator<Item = io::Result<Box<dyn Series + '_>>> + '_> {
        (**self).series()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = Box<dyn Sample + '_>> + '_> {
        (**self).iter()
    }
}
