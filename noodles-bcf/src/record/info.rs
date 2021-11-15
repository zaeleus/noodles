/// BCF record info.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Info {
    buf: Vec<u8>,
    field_count: usize,
}

impl Info {
    /// Returns the number of info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Info;
    /// let info = Info::default();
    /// assert_eq!(info.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.field_count
    }

    pub(crate) fn set_field_count(&mut self, field_count: usize) {
        self.field_count = field_count;
    }
}

impl AsRef<[u8]> for Info {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl AsMut<Vec<u8>> for Info {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }
}
