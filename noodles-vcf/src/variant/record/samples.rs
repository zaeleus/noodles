//! Variant record samples.

mod sample;
pub mod series;

pub use self::{sample::Sample, series::Series};

/// Variant record samples.
pub trait Samples {
    /// Returns whether there are any samples.
    fn is_empty(&self) -> bool;

    /// Returns an iterator over series.
    fn series(&self) -> Box<dyn Iterator<Item = Box<dyn Series + '_>> + '_>;

    /// Returns an iterator over samples.
    fn samples(&self) -> Box<dyn Iterator<Item = Box<dyn Sample + '_>> + '_>;
}
