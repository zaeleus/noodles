use std::collections::HashMap;

use crate::container::block;

const LOW_READER_COUNT: usize = 64;
const LOW_READER_COUNT_I32: i32 = LOW_READER_COUNT as i32;

pub struct ExternalDataReaders<'c> {
    low_readers: [Option<&'c [u8]>; LOW_READER_COUNT],
    high_readers: HashMap<block::ContentId, &'c [u8]>,
}

impl<'c> ExternalDataReaders<'c> {
    pub fn new() -> Self {
        Self {
            low_readers: [None; LOW_READER_COUNT],
            high_readers: HashMap::new(),
        }
    }

    pub fn insert(&mut self, id: block::ContentId, reader: &'c [u8]) {
        match id {
            i @ 0..LOW_READER_COUNT_I32 => {
                self.low_readers[i as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &block::ContentId) -> Option<&mut &'c [u8]> {
        match *id {
            i @ 0..LOW_READER_COUNT_I32 => self.low_readers[i as usize].as_mut(),
            _ => self.high_readers.get_mut(id),
        }
    }
}
