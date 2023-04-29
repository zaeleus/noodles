use noodles_sam::record::data::field::{Tag, Type};

use crate::container::block;

/// A CRAM data container compression header preservation map tag IDs dictionary key.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Key {
    tag: Tag,
    ty: Type,
}

impl Key {
    /// Creates a tag IDs dictionary key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::data_container::compression_header::preservation_map::tag_ids_dictionary::Key;
    /// use noodles_sam::record::data::field::{tag, Type};
    /// let key = Key::new(tag::ALIGNMENT_HIT_COUNT, Type::UInt8);
    /// ```
    pub fn new(tag: Tag, ty: Type) -> Self {
        Self { tag, ty }
    }

    /// Returns the tag.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::data_container::compression_header::preservation_map::tag_ids_dictionary::Key;
    /// use noodles_sam::record::data::field::{tag, Type};
    /// let key = Key::new(tag::ALIGNMENT_HIT_COUNT, Type::UInt8);
    /// assert_eq!(key.tag(), tag::ALIGNMENT_HIT_COUNT);
    /// ```
    pub fn tag(self) -> Tag {
        self.tag
    }

    /// Returns the type.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_cram::data_container::compression_header::preservation_map::tag_ids_dictionary::Key;
    /// use noodles_sam::record::data::field::{tag, Type};
    /// let key = Key::new(tag::ALIGNMENT_HIT_COUNT, Type::UInt8);
    /// assert_eq!(key.ty(), Type::UInt8);
    /// ```
    pub fn ty(self) -> Type {
        self.ty
    }
}

impl From<Key> for block::ContentId {
    fn from(key: Key) -> Self {
        let [l, r] = key.tag.as_ref();
        let ty = u8::from(key.ty);
        let id = i32::from(*l) << 16 | i32::from(*r) << 8 | i32::from(ty);
        Self::from(id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_key_for_block_content_id() {
        use noodles_sam::record::data::field::tag;

        let key = Key::new(tag::COMMENT, Type::String);
        assert_eq!(block::ContentId::from(key), block::ContentId::from(4411226));

        let key = Key::new(tag::ALIGNMENT_HIT_COUNT, Type::Int32);
        assert_eq!(block::ContentId::from(key), block::ContentId::from(5130345));
    }
}
