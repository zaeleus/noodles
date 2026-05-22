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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let readers = ExternalDataReaders::new();
        assert!(readers.low_readers.iter().all(|reader| reader.is_none()));
        assert!(readers.high_readers.is_empty());
    }

    #[test]
    fn test_insert() {
        let mut readers = ExternalDataReaders::new();

        readers.insert(0, b"z");
        readers.insert(0, b"a");
        readers.insert(LOW_READER_COUNT_I32 - 1, b"b");
        readers.insert(LOW_READER_COUNT_I32, b"c");

        assert_eq!(readers.low_readers[0], Some(&b"a"[..]));
        assert_eq!(readers.low_readers[63], Some(&b"b"[..]));
        assert_eq!(readers.high_readers.get(&64), Some(&&b"c"[..]));
    }

    #[test]
    fn test_get_mut() {
        let mut readers = ExternalDataReaders::new();

        readers.insert(0, b"a");
        readers.insert(LOW_READER_COUNT_I32 - 1, b"b");
        readers.insert(LOW_READER_COUNT_I32, b"c");

        assert_eq!(readers.get_mut(&0), Some(&mut &b"a"[..]));
        assert_eq!(readers.get_mut(&63), Some(&mut &b"b"[..]));
        assert_eq!(readers.get_mut(&64), Some(&mut &b"c"[..]));
    }
}
