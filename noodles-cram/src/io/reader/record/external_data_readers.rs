use std::collections::HashMap;

use bytes::Buf;

use crate::container::block;

pub struct ExternalDataReaders<B> {
    low_readers: [Option<B>; 64],
    high_readers: HashMap<block::ContentId, B>,
}

impl<B> ExternalDataReaders<B>
where
    B: Buf,
{
    pub fn new() -> Self {
        Self {
            low_readers: init_low_readers(),
            high_readers: HashMap::new(),
        }
    }

    pub fn insert(&mut self, id: block::ContentId, reader: B) {
        match id {
            i @ 0..=63 => {
                self.low_readers[i as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &block::ContentId) -> Option<&mut B> {
        match *id {
            i @ 0..=63 => self.low_readers[i as usize].as_mut(),
            _ => self.high_readers.get_mut(id),
        }
    }
}

fn init_low_readers<B>() -> [Option<B>; 64]
where
    B: Buf,
{
    [
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None,
    ]
}
