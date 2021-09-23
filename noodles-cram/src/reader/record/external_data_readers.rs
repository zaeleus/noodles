use std::{collections::HashMap, io::Read};

use crate::num::Itf8;

pub struct ExternalDataReaders<R> {
    low_readers: [Option<R>; 64],
    high_readers: HashMap<Itf8, R>,
}

impl<R> ExternalDataReaders<R>
where
    R: Read,
{
    pub fn new() -> Self {
        Self {
            low_readers: init_low_readers(),
            high_readers: HashMap::new(),
        }
    }

    pub fn insert(&mut self, id: Itf8, reader: R) {
        match id {
            0..=63 => {
                self.low_readers[id as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &Itf8) -> Option<&mut R> {
        match *id {
            0..=63 => self.low_readers[*id as usize].as_mut(),
            _ => self.high_readers.get_mut(id),
        }
    }
}

fn init_low_readers<R>() -> [Option<R>; 64]
where
    R: Read,
{
    [
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None,
    ]
}
