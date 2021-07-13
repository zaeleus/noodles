use noodles_bam::record::data::field::value::Type;

/// A CRAM record tag key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Key {
    tag: [u8; 2],
    ty: Type,
}

impl Key {
    /// Creates a CRAM record tag key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::field::value::Type;
    /// use noodles_cram::record::tag::Key;
    /// let key = Key::new([b'N', b'H'], Type::Int8);
    /// ```
    pub fn new(tag: [u8; 2], ty: Type) -> Self {
        Self { tag, ty }
    }

    /// Returns the tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::field::value::Type;
    /// use noodles_cram::record::tag::Key;
    /// let key = Key::new([b'N', b'H'], Type::Int8);
    /// assert_eq!(key.tag(), [b'N', b'H']);
    /// ```
    pub fn tag(self) -> [u8; 2] {
        self.tag
    }

    /// Returns the type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::data::field::value::Type;
    /// use noodles_cram::record::tag::Key;
    /// let key = Key::new([b'N', b'H'], Type::Int8);
    /// assert_eq!(key.ty(), Type::Int8);
    /// ```
    pub fn ty(self) -> Type {
        self.ty
    }

    pub(crate) fn id(self) -> i32 {
        let [l, r] = self.tag;
        let ty = char::from(self.ty) as u8;
        i32::from(l) << 16 | i32::from(r) << 8 | i32::from(ty)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_id() {
        let key = Key::new([b'C', b'O'], Type::String);
        assert_eq!(key.id(), 4411226);

        let key = Key::new([b'N', b'H'], Type::Int32);
        assert_eq!(key.id(), 5130345);
    }
}
