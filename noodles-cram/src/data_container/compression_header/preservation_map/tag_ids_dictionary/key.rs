use noodles_sam::{
    self as sam,
    record::data::field::{value::Type, Tag},
};

use crate::container::block;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Key {
    tag: Tag,
    ty: Type,
}

impl Key {
    pub fn new(tag: Tag, ty: Type) -> Self {
        Self { tag, ty }
    }

    pub fn tag(self) -> Tag {
        self.tag
    }

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

impl From<&sam::record::data::Field> for Key {
    fn from(field: &sam::record::data::Field) -> Self {
        Self::new(field.tag(), field.value().ty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_key_for_block_content_id() {
        let key = Key::new(Tag::Comment, Type::String);
        assert_eq!(block::ContentId::from(key), block::ContentId::from(4411226));

        let key = Key::new(Tag::AlignmentHitCount, Type::Int32);
        assert_eq!(block::ContentId::from(key), block::ContentId::from(5130345));
    }
}
