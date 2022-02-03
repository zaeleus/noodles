//! FASTX record

/* std use */

/* crate use */

/* project use */
use crate::reader::Reader;
use crate::record::Record;

/* mod declaration */

pub struct Records<'a, R>
where
    R: std::io::BufRead,
{
    inner: &'a mut Reader<R>,
}

impl<'a, R> Records<'a, R>
where
    R: std::io::BufRead,
{
    pub fn new(inner: &'a mut Reader<R>) -> Self {
        Self { inner }
    }
}

impl<'a, R> Iterator for Records<'a, R>
where
    R: std::io::BufRead,
{
    type Item = std::io::Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next_record()
    }
}
