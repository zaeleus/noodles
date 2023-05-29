//! SAM header read group platform.

use std::{borrow::Cow, error, fmt, str::FromStr};

/// A SAM header read group platform (`PL`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Platform {
    /// Capillary electrophoresis sequencing (`CAPILLARY`).
    Capillary,
    /// DNBseq sequencing (`DNBSEQ`).
    DnbSeq,
    /// Element Biosciences (`ELEMENT`).
    Element,
    /// Helicos sequencing (`HELICOS`).
    Helicos,
    /// Illumina sequencing (`ILLUMINA`).
    Illumina,
    /// Ion Torrent sequencing (`IONTORRENT`).
    IonTorrent,
    /// 454 Life Sciences sequencing (`LS454`).
    Ls454,
    /// Oxford Nanopore Technologies (ONT) sequencing (`ONT`).
    Ont,
    /// Pacific Biosciences (PacBio) sequencing (`PACBIO`).
    PacBio,
    /// Singular Genomics (`SINGULAR`).
    Singular,
    /// SOLiD sequencing (`SOLID`).
    Solid,
    /// Ultima Genomics (`ULTIMA`).
    Ultima,
}

impl AsRef<str> for Platform {
    fn as_ref(&self) -> &str {
        match self {
            Self::Capillary => "CAPILLARY",
            Self::DnbSeq => "DNBSEQ",
            Self::Element => "ELEMENT",
            Self::Helicos => "HELICOS",
            Self::Illumina => "ILLUMINA",
            Self::IonTorrent => "IONTORRENT",
            Self::Ls454 => "LS454",
            Self::Ont => "ONT",
            Self::PacBio => "PACBIO",
            Self::Singular => "SINGULAR",
            Self::Solid => "SOLID",
            Self::Ultima => "ULTIMA",
        }
    }
}

impl fmt::Display for Platform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_ref())
    }
}

/// An error returned when a raw SAM header read group platform fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input is empty.
    Empty,
    /// The input is invalid.
    Invalid,
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("empty input"),
            Self::Invalid => f.write_str("invalid input"),
        }
    }
}

impl FromStr for Platform {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }

        let (is_lowercase, is_uppercase) = s
            .chars()
            .filter(char::is_ascii_alphabetic)
            .map(|c| c.is_ascii_uppercase())
            .fold((true, true), |(l, u), is_ascii_uppercase| {
                (l && !is_ascii_uppercase, u && is_ascii_uppercase)
            });

        let t = if is_uppercase {
            Cow::from(s)
        } else if is_lowercase {
            Cow::from(s.to_uppercase())
        } else {
            return Err(ParseError::Invalid);
        };

        match t.as_ref() {
            "CAPILLARY" => Ok(Self::Capillary),
            "DNBSEQ" => Ok(Self::DnbSeq),
            "ELEMENT" => Ok(Self::Element),
            "HELICOS" => Ok(Self::Helicos),
            "ILLUMINA" => Ok(Self::Illumina),
            "IONTORRENT" => Ok(Self::IonTorrent),
            "LS454" => Ok(Self::Ls454),
            "ONT" => Ok(Self::Ont),
            "PACBIO" => Ok(Self::PacBio),
            "SINGULAR" => Ok(Self::Singular),
            "SOLID" => Ok(Self::Solid),
            "ULTIMA" => Ok(Self::Ultima),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Platform::Capillary.to_string(), "CAPILLARY");
        assert_eq!(Platform::DnbSeq.to_string(), "DNBSEQ");
        assert_eq!(Platform::Element.to_string(), "ELEMENT");
        assert_eq!(Platform::Helicos.to_string(), "HELICOS");
        assert_eq!(Platform::Illumina.to_string(), "ILLUMINA");
        assert_eq!(Platform::IonTorrent.to_string(), "IONTORRENT");
        assert_eq!(Platform::Ls454.to_string(), "LS454");
        assert_eq!(Platform::Ont.to_string(), "ONT");
        assert_eq!(Platform::PacBio.to_string(), "PACBIO");
        assert_eq!(Platform::Singular.to_string(), "SINGULAR");
        assert_eq!(Platform::Solid.to_string(), "SOLID");
        assert_eq!(Platform::Ultima.to_string(), "ULTIMA");
    }

    #[test]
    fn test_from_str() {
        assert_eq!("CAPILLARY".parse(), Ok(Platform::Capillary));
        assert_eq!("DNBSEQ".parse(), Ok(Platform::DnbSeq));
        assert_eq!("ELEMENT".parse(), Ok(Platform::Element));
        assert_eq!("HELICOS".parse(), Ok(Platform::Helicos));
        assert_eq!("ILLUMINA".parse(), Ok(Platform::Illumina));
        assert_eq!("IONTORRENT".parse(), Ok(Platform::IonTorrent));
        assert_eq!("LS454".parse(), Ok(Platform::Ls454));
        assert_eq!("ONT".parse(), Ok(Platform::Ont));
        assert_eq!("PACBIO".parse(), Ok(Platform::PacBio));
        assert_eq!("SINGULAR".parse(), Ok(Platform::Singular));
        assert_eq!("SOLID".parse(), Ok(Platform::Solid));
        assert_eq!("ULTIMA".parse(), Ok(Platform::Ultima));

        assert_eq!("illumina".parse(), Ok(Platform::Illumina));

        assert_eq!("".parse::<Platform>(), Err(ParseError::Empty));
        assert_eq!("Illumina".parse::<Platform>(), Err(ParseError::Invalid));
        assert_eq!("NOODLES".parse::<Platform>(), Err(ParseError::Invalid));
    }
}
