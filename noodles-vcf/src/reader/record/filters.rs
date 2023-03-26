use std::{io, mem};

use indexmap::IndexSet;

use crate::record::Filters;

pub(super) fn parse_filters(s: &str, filters: &mut Option<Filters>) -> io::Result<()> {
    const DELIMITER: char = ';';
    const PASS: &str = "PASS";

    if s.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "empty filters"));
    } else if s == PASS {
        *filters = Some(Filters::Pass);
        return Ok(());
    }

    let mut set = match mem::take(filters) {
        Some(Filters::Pass) | None => IndexSet::new(),
        Some(Filters::Fail(mut set)) => {
            set.clear();
            set
        }
    };

    for raw_filter in s.split(DELIMITER) {
        if !set.insert(raw_filter.into()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "duplicate filter",
            ));
        } else if !is_valid_filter(raw_filter) {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid filter"));
        }
    }

    *filters = Some(Filters::Fail(set));

    Ok(())
}

fn is_valid_filter(s: &str) -> bool {
    match s {
        "" | "0" => false,
        _ => s.chars().all(|c| !c.is_ascii_whitespace()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_filters() -> io::Result<()> {
        let mut filters = None;

        parse_filters("PASS", &mut filters)?;
        assert_eq!(filters, Some(Filters::Pass));

        parse_filters("q10", &mut filters)?;
        assert_eq!(
            filters,
            Some(Filters::Fail([String::from("q10")].into_iter().collect()))
        );

        parse_filters("q10;s50", &mut filters)?;
        assert_eq!(
            filters,
            Some(Filters::Fail(
                [String::from("q10"), String::from("s50")]
                    .into_iter()
                    .collect()
            ))
        );

        assert!(matches!(
            parse_filters("", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters("q10;q10", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters("0", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters("q 10", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters(";q10", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters("q10;;s50", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));
        assert!(matches!(
            parse_filters("q10;", &mut filters),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
