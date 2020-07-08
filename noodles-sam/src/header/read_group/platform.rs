use std::{error, fmt, str::FromStr};

/// A SAM header read group platform (`PL`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Platform {
    /// Capillary electrophoresis sequencing (`CAPILLARY`).
    Capillary,
    /// DNBseq sequencing (`DNBSEQ`).
    DnbSeq,
    /// 454 Life Sciences sequencing (`LS454`).
    LS454,
    /// Illumina sequencing (`ILLUMINA`).
    Illumina,
    /// SOLiD sequencing (`SOLID`).
    Solid,
    /// Helicos sequencing (`HELICOS`).
    Helicos,
    /// Ion Torrent sequencing (`IONTORRENT`).
    IonTorrent,
    /// Oxford Nanopore Technologies (ONT) sequencing (`ONT`).
    Ont,
    /// PacBio sequencing (`PACBIO`).
    PacBio,
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
        match s {
            "" => Err(ParseError::Empty),
            "CAPILLARY" => Ok(Self::Capillary),
            "DNBSEQ" => Ok(Self::DnbSeq),
            "LS454" => Ok(Self::LS454),
            "ILLUMINA" => Ok(Self::Illumina),
            "SOLID" => Ok(Self::Solid),
            "HELICOS" => Ok(Self::Helicos),
            "IONTORRENT" => Ok(Self::IonTorrent),
            "ONT" => Ok(Self::Ont),
            "PACBIO" => Ok(Self::PacBio),
            _ => Err(ParseError::Invalid),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() {
        assert_eq!("CAPILLARY".parse(), Ok(Platform::Capillary));
        assert_eq!("DNBSEQ".parse(), Ok(Platform::DnbSeq));
        assert_eq!("LS454".parse(), Ok(Platform::LS454));
        assert_eq!("ILLUMINA".parse(), Ok(Platform::Illumina));
        assert_eq!("SOLID".parse(), Ok(Platform::Solid));
        assert_eq!("HELICOS".parse(), Ok(Platform::Helicos));
        assert_eq!("IONTORRENT".parse(), Ok(Platform::IonTorrent));
        assert_eq!("ONT".parse(), Ok(Platform::Ont));
        assert_eq!("PACBIO".parse(), Ok(Platform::PacBio));

        assert_eq!("".parse::<Platform>(), Err(ParseError::Empty));
        assert_eq!("NOODLES".parse::<Platform>(), Err(ParseError::Invalid));
        assert_eq!("Illumina".parse::<Platform>(), Err(ParseError::Invalid));
        assert_eq!("illumina".parse::<Platform>(), Err(ParseError::Invalid));
    }
}
