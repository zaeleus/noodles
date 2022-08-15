use std::io::{self, Read};

use crate::Block;

pub struct Reader<R> {
    inner: R,
    buf: Vec<u8>,
}

impl<R> Reader<R>
where
    R: Read,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
        }
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
        use super::{parse_frame, read_frame_into};

        if read_frame_into(&mut self.inner, &mut self.buf)?.is_some() {
            parse_frame(&self.buf).map(Some)
        } else {
            Ok(None)
        }
    }
}
