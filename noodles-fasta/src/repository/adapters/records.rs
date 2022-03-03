use std::io;

use crate::{repository::Adapter, Record};

impl Adapter for Vec<Record> {
    fn get(&mut self, name: &str) -> Option<io::Result<Record>> {
        self.iter()
            .find(|record| record.name() == name)
            .cloned()
            .map(Ok)
    }
}
