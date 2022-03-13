use noodles_sam::record::data::field::{value::Type, Tag};

/// A CRAM record tag key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Key {
    tag: Tag,
    ty: Type,
}

impl Key {
    /// Creates a CRAM record tag key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::tag::Key;
    /// use noodles_sam::record::data::field::{value::Type, Tag};
    /// let key = Key::new(Tag::AlignmentHitCount, Type::Int8);
    /// ```
    pub fn new(tag: Tag, ty: Type) -> Self {
        Self { tag, ty }
    }

    /// Returns the tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::tag::Key;
    /// use noodles_sam::record::data::field::{value::Type, Tag};
    /// let key = Key::new(Tag::AlignmentHitCount, Type::Int8);
    /// assert_eq!(key.tag(), Tag::AlignmentHitCount);
    /// ```
    pub fn tag(self) -> Tag {
        self.tag
    }

    /// Returns the type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::record::tag::Key;
    /// use noodles_sam::record::data::field::{value::Type, Tag};
    /// let key = Key::new(Tag::AlignmentHitCount, Type::Int8);
    /// assert_eq!(key.ty(), Type::Int8);
    /// ```
    pub fn ty(self) -> Type {
        self.ty
    }

    pub(crate) fn id(self) -> i32 {
        let [l, r] = self.tag.as_ref();
        let ty = u8::from(self.ty);
        i32::from(*l) << 16 | i32::from(*r) << 8 | i32::from(ty)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_id() {
        let key = Key::new(Tag::Comment, Type::String);
        assert_eq!(key.id(), 4411226);

        let key = Key::new(Tag::AlignmentHitCount, Type::Int32);
        assert_eq!(key.id(), 5130345);
    }
}
