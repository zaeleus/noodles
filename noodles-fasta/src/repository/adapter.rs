use std::io;

use crate::Record;

/// A repository adapter.
pub trait Adapter: Send + Sync {
    /// Returns the record with the given name.
    fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>>;
}
