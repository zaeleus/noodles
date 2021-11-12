#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Genotypes {
    buf: Vec<u8>,
    format_count: usize,
    sample_count: usize,
}

impl Genotypes {
    pub(crate) fn len(&self) -> usize {
        self.sample_count
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(crate) fn format_count(&self) -> usize {
        self.format_count
    }

    pub(crate) fn set_format_count(&mut self, format_count: usize) {
        self.format_count = format_count;
    }

    pub(crate) fn set_sample_count(&mut self, sample_count: usize) {
        self.sample_count = sample_count;
    }
}

impl AsRef<[u8]> for Genotypes {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl AsMut<Vec<u8>> for Genotypes {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }
}
