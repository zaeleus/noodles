use std::fmt;

use noodles_sam as sam;

/// A BAM record sequence base.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    /// Equal to the reference base.
    Eq,
    /// Adenine.
    A,
    /// Cytosine.
    C,
    /// Amino.
    M,
    /// Guanine.
    G,
    /// Purine.
    R,
    /// Strong.
    S,
    /// Not T.
    V,
    /// Thymine.
    T,
    /// Weak.
    W,
    /// Pyrimidine.
    Y,
    /// Not G.
    H,
    /// Keto.
    K,
    /// Not C.
    D,
    /// Not A.
    B,
    /// Any base.
    N,
}

impl Base {
    /// Returns the complement of this base.
    ///
    /// See <https://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation>.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bam::record::sequence::Base;
    /// assert_eq!(Base::A.complement(), Base::T);
    /// assert_eq!(Base::C.complement(), Base::G);
    /// assert_eq!(Base::G.complement(), Base::C);
    /// assert_eq!(Base::T.complement(), Base::A);
    /// ```
    pub fn complement(self) -> Self {
        match self {
            Self::Eq => Self::Eq,
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::W => Self::W,
            Self::S => Self::S,
            Self::M => Self::K,
            Self::K => Self::M,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::B => Self::V,
            Self::D => Self::H,
            Self::H => Self::D,
            Self::V => Self::B,
            Self::N => Self::N,
        }
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

impl From<sam::record::sequence::Base> for Base {
    fn from(base: sam::record::sequence::Base) -> Self {
        use sam::record::sequence::Base as SamBase;

        match base {
            SamBase::Eq => Self::Eq,
            SamBase::A => Self::A,
            SamBase::C => Self::C,
            SamBase::M => Self::M,
            SamBase::G => Self::G,
            SamBase::R => Self::R,
            SamBase::S => Self::S,
            SamBase::V => Self::V,
            SamBase::T => Self::T,
            SamBase::W => Self::W,
            SamBase::Y => Self::Y,
            SamBase::H => Self::H,
            SamBase::K => Self::K,
            SamBase::D => Self::D,
            SamBase::B => Self::B,
            // § 4.2.3 SEQ and QUAL encoding (2020-06-19): "all other characters [map] to 'N'".
            _ => Self::N,
        }
    }
}

impl From<Base> for char {
    fn from(base: Base) -> Self {
        match base {
            Base::Eq => '=',
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::W => 'W',
            Base::S => 'S',
            Base::M => 'M',
            Base::K => 'K',
            Base::R => 'R',
            Base::Y => 'Y',
            Base::B => 'B',
            Base::D => 'D',
            Base::H => 'H',
            Base::V => 'V',
            Base::N => 'N',
        }
    }
}

impl From<Base> for u8 {
    fn from(base: Base) -> Self {
        match base {
            Base::Eq => 0,
            Base::A => 1,
            Base::C => 2,
            Base::M => 3,
            Base::G => 4,
            Base::R => 5,
            Base::S => 6,
            Base::V => 7,
            Base::T => 8,
            Base::W => 9,
            Base::Y => 10,
            Base::H => 11,
            Base::K => 12,
            Base::D => 13,
            Base::B => 14,
            Base::N => 15,
        }
    }
}

impl From<Base> for sam::record::sequence::Base {
    fn from(base: Base) -> Self {
        match base {
            Base::Eq => Self::Eq,
            Base::A => Self::A,
            Base::C => Self::C,
            Base::G => Self::G,
            Base::T => Self::T,
            Base::W => Self::W,
            Base::S => Self::S,
            Base::M => Self::M,
            Base::K => Self::K,
            Base::R => Self::R,
            Base::Y => Self::Y,
            Base::B => Self::B,
            Base::D => Self::D,
            Base::H => Self::H,
            Base::V => Self::V,
            Base::N => Self::N,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(Base::Eq.complement(), Base::Eq);
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::M.complement(), Base::K);
        assert_eq!(Base::G.complement(), Base::C);
        assert_eq!(Base::R.complement(), Base::Y);
        assert_eq!(Base::S.complement(), Base::S);
        assert_eq!(Base::V.complement(), Base::B);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::W.complement(), Base::W);
        assert_eq!(Base::Y.complement(), Base::R);
        assert_eq!(Base::H.complement(), Base::D);
        assert_eq!(Base::K.complement(), Base::M);
        assert_eq!(Base::D.complement(), Base::H);
        assert_eq!(Base::B.complement(), Base::V);
        assert_eq!(Base::N.complement(), Base::N);
    }

    #[test]
    fn test_fmt() {
        assert_eq!(Base::Eq.to_string(), "=");
        assert_eq!(Base::A.to_string(), "A");
        assert_eq!(Base::C.to_string(), "C");
        assert_eq!(Base::M.to_string(), "M");
        assert_eq!(Base::G.to_string(), "G");
        assert_eq!(Base::R.to_string(), "R");
        assert_eq!(Base::S.to_string(), "S");
        assert_eq!(Base::V.to_string(), "V");
        assert_eq!(Base::T.to_string(), "T");
        assert_eq!(Base::W.to_string(), "W");
        assert_eq!(Base::Y.to_string(), "Y");
        assert_eq!(Base::H.to_string(), "H");
        assert_eq!(Base::K.to_string(), "K");
        assert_eq!(Base::D.to_string(), "D");
        assert_eq!(Base::B.to_string(), "B");
        assert_eq!(Base::N.to_string(), "N");
    }

    #[test]
    fn test_from_sam_record_sequence_base_for_base() {
        use sam::record::sequence::Base as SamBase;

        assert_eq!(Base::from(SamBase::Eq), Base::Eq);
        assert_eq!(Base::from(SamBase::A), Base::A);
        assert_eq!(Base::from(SamBase::C), Base::C);
        assert_eq!(Base::from(SamBase::G), Base::G);
        assert_eq!(Base::from(SamBase::T), Base::T);
        assert_eq!(Base::from(SamBase::W), Base::W);
        assert_eq!(Base::from(SamBase::S), Base::S);
        assert_eq!(Base::from(SamBase::M), Base::M);
        assert_eq!(Base::from(SamBase::K), Base::K);
        assert_eq!(Base::from(SamBase::R), Base::R);
        assert_eq!(Base::from(SamBase::Y), Base::Y);
        assert_eq!(Base::from(SamBase::B), Base::B);
        assert_eq!(Base::from(SamBase::D), Base::D);
        assert_eq!(Base::from(SamBase::H), Base::H);
        assert_eq!(Base::from(SamBase::V), Base::V);
        assert_eq!(Base::from(SamBase::N), Base::N);

        assert_eq!(Base::from(SamBase::Z), Base::N);
    }

    #[test]
    fn test_from_base_for_char() {
        assert_eq!(char::from(Base::Eq), '=');
        assert_eq!(char::from(Base::A), 'A');
        assert_eq!(char::from(Base::C), 'C');
        assert_eq!(char::from(Base::M), 'M');
        assert_eq!(char::from(Base::G), 'G');
        assert_eq!(char::from(Base::R), 'R');
        assert_eq!(char::from(Base::S), 'S');
        assert_eq!(char::from(Base::V), 'V');
        assert_eq!(char::from(Base::T), 'T');
        assert_eq!(char::from(Base::W), 'W');
        assert_eq!(char::from(Base::Y), 'Y');
        assert_eq!(char::from(Base::H), 'H');
        assert_eq!(char::from(Base::K), 'K');
        assert_eq!(char::from(Base::D), 'D');
        assert_eq!(char::from(Base::B), 'B');
        assert_eq!(char::from(Base::N), 'N');
    }

    #[test]
    fn test_from_base_for_u8() {
        assert_eq!(u8::from(Base::Eq), 0);
        assert_eq!(u8::from(Base::A), 1);
        assert_eq!(u8::from(Base::C), 2);
        assert_eq!(u8::from(Base::M), 3);
        assert_eq!(u8::from(Base::G), 4);
        assert_eq!(u8::from(Base::R), 5);
        assert_eq!(u8::from(Base::S), 6);
        assert_eq!(u8::from(Base::V), 7);
        assert_eq!(u8::from(Base::T), 8);
        assert_eq!(u8::from(Base::W), 9);
        assert_eq!(u8::from(Base::Y), 10);
        assert_eq!(u8::from(Base::H), 11);
        assert_eq!(u8::from(Base::K), 12);
        assert_eq!(u8::from(Base::D), 13);
        assert_eq!(u8::from(Base::B), 14);
        assert_eq!(u8::from(Base::N), 15);
    }

    #[test]
    fn test_from_base_for_sam_record_sequence_base() {
        use sam::record::sequence::Base as SamBase;

        assert_eq!(SamBase::from(Base::Eq), SamBase::Eq);
        assert_eq!(SamBase::from(Base::A), SamBase::A);
        assert_eq!(SamBase::from(Base::C), SamBase::C);
        assert_eq!(SamBase::from(Base::M), SamBase::M);
        assert_eq!(SamBase::from(Base::G), SamBase::G);
        assert_eq!(SamBase::from(Base::R), SamBase::R);
        assert_eq!(SamBase::from(Base::S), SamBase::S);
        assert_eq!(SamBase::from(Base::V), SamBase::V);
        assert_eq!(SamBase::from(Base::T), SamBase::T);
        assert_eq!(SamBase::from(Base::W), SamBase::W);
        assert_eq!(SamBase::from(Base::Y), SamBase::Y);
        assert_eq!(SamBase::from(Base::H), SamBase::H);
        assert_eq!(SamBase::from(Base::K), SamBase::K);
        assert_eq!(SamBase::from(Base::D), SamBase::D);
        assert_eq!(SamBase::from(Base::B), SamBase::B);
        assert_eq!(SamBase::from(Base::N), SamBase::N);
    }
}
