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
    pub fn get(&self, reference_base: Base, substitution_code: u8) -> Base {
        let i = match reference_base {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
            Base::N => 4,
        };

        let j = usize::from(substitution_code & 0x03);

        self.0[i][j]
    }

    pub fn find_code(&self, reference_base: Base, read_base: Base) -> u8 {
        for code in [0b00, 0b01, 0b10, 0b11] {
            if self.get(reference_base, code) == read_base {
                return code;
            }
        }

        unreachable!();
    }
}

impl Default for SubstitutionMatrix {
    fn default() -> Self {
        SubstitutionMatrix([
            [Base::C, Base::G, Base::T, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::A, Base::C, Base::T, Base::N],
            [Base::A, Base::C, Base::G, Base::N],
            [Base::A, Base::C, Base::G, Base::T],
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_code() {
        let matrix = SubstitutionMatrix([
            [Base::T, Base::G, Base::C, Base::N],
            [Base::A, Base::G, Base::T, Base::N],
            [Base::N, Base::A, Base::C, Base::T],
            [Base::G, Base::N, Base::A, Base::C],
            [Base::C, Base::G, Base::T, Base::A],
        ]);

        assert_eq!(matrix.find_code(Base::A, Base::T), 0b00);
        assert_eq!(matrix.find_code(Base::C, Base::G), 0b01);
        assert_eq!(matrix.find_code(Base::G, Base::C), 0b10);
        assert_eq!(matrix.find_code(Base::T, Base::C), 0b11);
    }
}
