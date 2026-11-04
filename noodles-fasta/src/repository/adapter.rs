use std::io;

use noodles_core::region::Interval;

use crate::{Record, record::Sequence};

/// A repository adapter.
pub trait Adapter: Send + Sync {
    /// Returns the record with the given name.
    fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>>;

    /// Returns the subsequence in the given interval of the sequence with the given name.
    ///
    /// This is optional, for adapters that can read part of a sequence without reading all of
    /// it. The subsequence ends early if the interval extends past the end of the sequence.
    ///
    /// The default implementation returns `None`, i.e., ranges are unsupported, and callers use
    /// [`Self::get`] instead.
    fn get_range(&mut self, name: &[u8], interval: Interval) -> Option<io::Result<Sequence>> {
        let _ = (name, interval);
        None
    }
}
