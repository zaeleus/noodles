use std::str::FromStr;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Platform {
    Capillary,
    LS454,
    Illumina,
    Solid,
    Helicos,
    IonTorrent,
    Ont,
    PacBio,
}

impl FromStr for Platform {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "CAPILLARY" => Ok(Self::Capillary),
            "LS454" => Ok(Self::LS454),
            "ILLUMINA" => Ok(Self::Illumina),
            "SOLID" => Ok(Self::Solid),
            "HELICOS" => Ok(Self::Helicos),
            "IONTORRENT" => Ok(Self::IonTorrent),
            "ONT" => Ok(Self::Ont),
            "PACBIO" => Ok(Self::PacBio),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ()> {
        assert_eq!("CAPILLARY".parse::<Platform>()?, Platform::Capillary);
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
