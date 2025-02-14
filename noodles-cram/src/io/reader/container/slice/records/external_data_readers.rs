use std::collections::HashMap;

use crate::container::block;

pub struct ExternalDataReaders<'c> {
    low_readers: [Option<&'c [u8]>; 64],
    high_readers: HashMap<block::ContentId, &'c [u8]>,
}

impl<'c> ExternalDataReaders<'c> {
    pub fn new() -> Self {
        Self {
            low_readers: init_low_readers(),
            high_readers: HashMap::new(),
        }
    }

    pub fn insert(&mut self, id: block::ContentId, reader: &'c [u8]) {
        match id {
            i @ 0..=63 => {
                self.low_readers[i as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &block::ContentId) -> Option<&mut &'c [u8]> {
        match *id {
            i @ 0..=63 => self.low_readers[i as usize].as_mut(),
            _ => self.high_readers.get_mut(id),
        }
    }
}

fn init_low_readers<'c>() -> [Option<&'c [u8]>; 64] {
    [
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None,
    ]
}
