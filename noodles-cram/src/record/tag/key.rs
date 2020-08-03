use noodles_bam::record::data::field::value::Type;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Key {
    tag: [u8; 2],
    ty: Type,
}

impl Key {
    pub fn new(tag: [u8; 2], ty: Type) -> Self {
        Self { tag, ty }
    }

    pub fn tag(&self) -> [u8; 2] {
        self.tag
    }

    pub fn ty(&self) -> Type {
        self.ty
    }

    pub fn id(&self) -> i32 {
        let [l, r] = self.tag;
        let ty = char::from(self.ty) as u8;
        (l as i32) << 16 | (r as i32) << 8 | (ty as i32)
    }
}
