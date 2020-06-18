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
#[derive(Debug, Eq, PartialEq)]
pub struct ParseError(String);

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid platform: {}", self.0)
    }
}

impl FromStr for Platform {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "CAPILLARY" => Ok(Self::Capillary),
            "DNBSEQ" => Ok(Self::DnbSeq),
            "LS454" => Ok(Self::LS454),
            "ILLUMINA" => Ok(Self::Illumina),
            "SOLID" => Ok(Self::Solid),
            "HELICOS" => Ok(Self::Helicos),
            "IONTORRENT" => Ok(Self::IonTorrent),
            "ONT" => Ok(Self::Ont),
            "PACBIO" => Ok(Self::PacBio),
            _ => Err(ParseError(s.into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        assert_eq!("CAPILLARY".parse::<Platform>()?, Platform::Capillary);
        assert_eq!("DNBSEQ".parse::<Platform>()?, Platform::DnbSeq);
        assert_eq!("LS454".parse::<Platform>()?, Platform::LS454);
        assert_eq!("ILLUMINA".parse::<Platform>()?, Platform::Illumina);
        assert_eq!("SOLID".parse::<Platform>()?, Platform::Solid);
        assert_eq!("HELICOS".parse::<Platform>()?, Platform::Helicos);
        assert_eq!("IONTORRENT".parse::<Platform>()?, Platform::IonTorrent);
        assert_eq!("ONT".parse::<Platform>()?, Platform::Ont);
        assert_eq!("PACBIO".parse::<Platform>()?, Platform::PacBio);

        assert!("".parse::<Platform>().is_err());
        assert!("NOODLES".parse::<Platform>().is_err());
        assert!("Illumina".parse::<Platform>().is_err());
        assert!("illumina".parse::<Platform>().is_err());

        Ok(())
    }
}
