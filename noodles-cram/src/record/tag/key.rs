use noodles_bam::record::data::field::value::Type;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Key {
    tag: [u8; 2],
    ty: Type,
}

impl Key {
    pub fn new(tag: [u8; 2], ty: Type) -> Self {
        Self { tag, ty }
    }

    pub fn tag(self) -> [u8; 2] {
        self.tag
    }

    pub fn ty(self) -> Type {
        self.ty
    }

    pub fn id(self) -> i32 {
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
