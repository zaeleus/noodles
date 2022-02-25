use std::io;

use crate::Record;

/// A repository adapter.
pub trait Adapter {
    /// Returns the record with the given name.
    fn get(&mut self, name: &str) -> Option<io::Result<Record>>;
}

impl Adapter for Vec<Record> {
    fn get(&mut self, name: &str) -> Option<io::Result<Record>> {
        self.iter()
            .find(|record| record.name() == name)
            .cloned()
            .map(Ok)
    }
}
