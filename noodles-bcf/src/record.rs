//! BCF record and fields.

use std::ops::{Deref, DerefMut};

/// A BCF record.
#[derive(Clone, Default, Eq, PartialEq)]
pub struct Record(Vec<u8>);

impl Record {
    pub(crate) fn resize(&mut self, new_len: usize) {
        self.0.resize(new_len, Default::default());
    }
}

impl Deref for Record {
    type Target = [u8];

    fn deref(&self) -> &[u8] {
        &self.0
    }
}

impl DerefMut for Record {
    fn deref_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}
