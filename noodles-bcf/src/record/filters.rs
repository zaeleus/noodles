#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Filters(Vec<usize>);

impl Filters {
    pub(crate) fn len(&self) -> usize {
        self.0.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl AsRef<[usize]> for Filters {
    fn as_ref(&self) -> &[usize] {
        &self.0
    }
}

impl AsMut<Vec<usize>> for Filters {
    fn as_mut(&mut self) -> &mut Vec<usize> {
        &mut self.0
    }
}
