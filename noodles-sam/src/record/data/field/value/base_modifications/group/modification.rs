use std::{error, fmt};

/// 5-Methylcytosine (m; 5mC; CHEBI:27551).
const FIVE_METHYLCYTOSINE: Modification = Modification::FiveMethylcytosine;

/// 5-Hydroxymethylcytosine (h; 5hmC; CHEBI:76792).
const FIVE_HYDROXYMETHYLCYTOSINE: Modification = Modification::FiveHydroxymethylcytosine;

/// 5-Formylcytosine (f; 5fC; CHEBI:76794).
const FIVE_FORMYLCYTOSINE: Modification = Modification::FiveFormylcytosine;

/// 5-Carboxylcytosine (c; 5caC; CHEBI:76793).
const FIVE_CARBOXYLCYTOSINE: Modification = Modification::FiveCarboxylcytosine;

/// 5-Hydroxymethyluracil (g; 5hmU; CHEBI:16964).
const FIVE_HYDROXYMETHYLURACIL: Modification = Modification::FiveHydroxymethyluracil;

/// 5-Formyluracil (e; 5fU; CHEBI:80961).
const FIVE_FORMYLURACIL: Modification = Modification::FiveFormyluracil;

/// 5-Carboxyluracil (b; 5caU; CHEBI:17477).
const FIVE_CARBOXYLURACIL: Modification = Modification::FiveCarboxyluracil;

/// 6-Methyladenine (a; 6mA; CHEBI:28871).
const SIX_METHYLADENINE: Modification = Modification::SixMethyladenine;

/// 8-Oxoguanine (o; 8oxoG; CHEBI:44605).
const EIGHT_OXOGUANINE: Modification = Modification::EightOxoguanine;

/// Xanthosine (n; Xao; CHEBI:18107).
const XANTHOSINE: Modification = Modification::Xanthosine;

/// A base modification.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Modification {
    /// 5-Methylcytosine (m; 5mC; CHEBI:27551)
    FiveMethylcytosine,
    /// 5-Hydroxymethylcytosine (h; 5hmC; CHEBI:76792)
    FiveHydroxymethylcytosine,
    /// 5-Formylcytosine (f; 5fC; CHEBI:76794)
    FiveFormylcytosine,
    /// 5-Carboxylcytosine (c; 5caC; CHEBI:76793)
    FiveCarboxylcytosine,
    /// 5-Hydroxymethyluracil (g; 5hmU; CHEBI:16964)
    FiveHydroxymethyluracil,
    /// 5-Formyluracil (e; 5fU; CHEBI:80961)
    FiveFormyluracil,
    /// 5-Carboxyluracil (b; 5caU; CHEBI:17477)
    FiveCarboxyluracil,
    /// 6-Methyladenine (a; 6mA; CHEBI:28871)
    SixMethyladenine,
    /// 8-Oxoguanine (o; 8oxoG; CHEBI:44605)
    EightOxoguanine,
    // Xanthosine (n; Xao; CHEBI:18107)
    Xanthosine,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid => write!(f, "invalid input"),
        }
    }
}

impl TryFrom<u8> for Modification {
    type Error = ParseError;

    fn try_from(b: u8) -> Result<Self, Self::Error> {
        match b {
            b'm' => Ok(Self::FiveMethylcytosine),
            b'h' => Ok(Self::FiveHydroxymethylcytosine),
            b'f' => Ok(Self::FiveFormylcytosine),
            b'c' => Ok(Self::FiveCarboxylcytosine),
            b'g' => Ok(Self::FiveHydroxymethyluracil),
            b'e' => Ok(Self::FiveFormyluracil),
            b'b' => Ok(Self::FiveCarboxyluracil),
            b'a' => Ok(Self::SixMethyladenine),
            b'o' => Ok(Self::EightOxoguanine),
            b'n' => Ok(Self::Xanthosine),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_u8_for_modification() {
        fn t(b: u8, expected: Modification) {
            assert_eq!(Modification::try_from(b), Ok(expected));
        }

        t(b'm', Modification::FiveMethylcytosine);
        t(b'h', Modification::FiveHydroxymethylcytosine);
        t(b'f', Modification::FiveFormylcytosine);
        t(b'c', Modification::FiveCarboxylcytosine);
        t(b'g', Modification::FiveHydroxymethyluracil);
        t(b'e', Modification::FiveFormyluracil);
        t(b'b', Modification::FiveCarboxyluracil);
        t(b'a', Modification::SixMethyladenine);
        t(b'o', Modification::EightOxoguanine);
        t(b'n', Modification::Xanthosine);

        assert_eq!(Modification::try_from(b'?'), Err(ParseError::Invalid));
    }
}
