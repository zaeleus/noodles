use std::io;

use crate::{Record, repository::Adapter};

impl Adapter for Vec<Record> {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>> {
        self.iter()
            .find(|record| record.name() == name)
            .cloned()
            .map(Ok)
    }
}
