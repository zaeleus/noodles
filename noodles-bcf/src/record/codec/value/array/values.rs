use std::io;

pub trait Values<'a, N> {
    fn len(&self) -> usize;
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_>;
}

impl<N> Values<'_, N> for Vec<N>
where
    N: Copy,
{
    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        Box::new((**self).iter().copied().map(Ok))
    }
}
