use std::str::FromStr;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GroupOrder {
    None,
    Query,
    Reference,
}

impl Default for GroupOrder {
    fn default() -> Self {
        Self::None
    }
}

impl FromStr for GroupOrder {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "query" => Ok(Self::Query),
            "reference" => Ok(Self::Reference),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(GroupOrder::default(), GroupOrder::None);
    }

    #[test]
    fn test_from_str() -> Result<(), ()> {
        assert_eq!("none".parse::<GroupOrder>()?, GroupOrder::None);
        assert_eq!("query".parse::<GroupOrder>()?, GroupOrder::Query);
        assert_eq!("reference".parse::<GroupOrder>()?, GroupOrder::Reference);

        assert!("".parse::<GroupOrder>().is_err());
        assert!("noodles".parse::<GroupOrder>().is_err());
        assert!("Query".parse::<GroupOrder>().is_err());

        Ok(())
    }
}
