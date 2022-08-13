use std::io::{self, Read};

use crate::Block;

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    pub fn get_ref(&self) -> &R {
        &self.inner
    }

    pub fn get_mut(&mut self) -> &mut R {
        &mut self.inner
    }

    pub fn into_inner(self) -> R {
        self.inner
    }

    pub fn next_block(&mut self) -> io::Result<Option<Block>> {
        use super::{parse_frame, read_frame};

        match read_frame(&mut self.inner)? {
            Some(src) => parse_frame(&src).map(Some),
            None => Ok(None),
        }
    }
}
