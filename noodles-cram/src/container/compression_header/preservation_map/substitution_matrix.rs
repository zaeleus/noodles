mod base;

pub use self::base::Base;

pub const READ_BASES: [[Base; 4]; 5] = [
    [Base::C, Base::G, Base::T, Base::N], // A
    [Base::A, Base::G, Base::T, Base::N], // C
    [Base::A, Base::C, Base::T, Base::N], // G
    [Base::A, Base::C, Base::G, Base::N], // T
    [Base::A, Base::C, Base::G, Base::T], // N
];

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SubstitutionMatrix(pub(crate) [[Base; 4]; 5]);

impl SubstitutionMatrix {
    pub fn get(&self, reference_base: Base, code: u8) -> Base {
        let i = reference_base as usize;
        let j = usize::from(code & 0x03);
        self.0[i][j]
    }

    pub fn find(&self, reference_base: Base, read_base: Base) -> u8 {
        const CODES: [u8; 4] = [0b00, 0b01, 0b10, 0b11];

        CODES
            .into_iter()
            .find(|code| self.get(reference_base, *code) == read_base)
            .unwrap()
    }
}

impl Default for SubstitutionMatrix {
    fn default() -> Self {
        SubstitutionMatrix(READ_BASES)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find() {
        let matrix = SubstitutionMatrix([
            [Base::T, Base::G, Base::C, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::N, Base::A, Base::C, Base::T],
            [Base::G, Base::N, Base::A, Base::C],
            [Base::C, Base::G, Base::T, Base::A],
        ]);

        assert_eq!(matrix.find(Base::A, Base::T), 0b00);
        assert_eq!(matrix.find(Base::C, Base::G), 0b01);
        assert_eq!(matrix.find(Base::G, Base::C), 0b10);
        assert_eq!(matrix.find(Base::T, Base::C), 0b11);
    }
}
