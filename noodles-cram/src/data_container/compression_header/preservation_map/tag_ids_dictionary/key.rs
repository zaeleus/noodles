use noodles_sam::{
    self as sam,
    record::data::field::{value::Type, Tag},
};

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

    pub fn id(self) -> i32 {
        let [l, r] = self.tag.as_ref();
        let ty = u8::from(self.ty);
        i32::from(*l) << 16 | i32::from(*r) << 8 | i32::from(ty)
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
    fn test_id() {
        let key = Key::new(Tag::Comment, Type::String);
        assert_eq!(key.id(), 4411226);

        let key = Key::new(Tag::AlignmentHitCount, Type::Int32);
        assert_eq!(key.id(), 5130345);
    }
}
