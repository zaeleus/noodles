use std::collections::HashMap;

use bytes::Buf;

pub struct ExternalDataReaders<B> {
    low_readers: [Option<B>; 64],
    high_readers: HashMap<i32, B>,
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

    pub fn insert(&mut self, id: i32, reader: B) {
        match id {
            0..=63 => {
                self.low_readers[id as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &i32) -> Option<&mut B> {
        match *id {
            0..=63 => self.low_readers[*id as usize].as_mut(),
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
