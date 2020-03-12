use std::str::FromStr;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SubsortOrder {
    Unsorted(String),
    QueryName(String),
    Coordinate(String),
}

impl FromStr for SubsortOrder {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut pieces = s.splitn(2, ':');

        let order = pieces.next().ok_or_else(|| ())?;
        let subsort = pieces.next().map(|s| s.into()).ok_or_else(|| ())?;

        match order {
            "unsorted" => Ok(Self::Unsorted(subsort)),
            "queryname" => Ok(Self::QueryName(subsort)),
            "coordinate" => Ok(Self::Coordinate(subsort)),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ()> {
        assert_eq!(
            "unsorted:MI".parse::<SubsortOrder>()?,
            SubsortOrder::Unsorted(String::from("MI"))
        );

        assert_eq!(
            "queryname:MI".parse::<SubsortOrder>()?,
            SubsortOrder::QueryName(String::from("MI"))
        );

        assert_eq!(
            "coordinate:MI".parse::<SubsortOrder>()?,
            SubsortOrder::Coordinate(String::from("MI"))
        );

        assert_eq!(
            "unsorted:MI:coordinate".parse::<SubsortOrder>()?,
            SubsortOrder::Unsorted(String::from("MI:coordinate"))
        );

        assert!("".parse::<SubsortOrder>().is_err());
        assert!("noodles".parse::<SubsortOrder>().is_err());
        assert!("queryname".parse::<SubsortOrder>().is_err());
        assert!("QueryName".parse::<SubsortOrder>().is_err());

        Ok(())
    }
}
