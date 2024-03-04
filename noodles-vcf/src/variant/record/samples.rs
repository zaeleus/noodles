//! Variant record samples.

mod sample;
pub mod series;

use std::io;

pub use self::{sample::Sample, series::Series};
use crate::Header;

/// Variant record samples.
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

    fn series(&self) -> Box<dyn Iterator<Item = io::Result<Box<dyn Series + '_>>> + '_> {
        (**self).series()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = Box<dyn Sample + '_>> + '_> {
        (**self).iter()
    }
}
