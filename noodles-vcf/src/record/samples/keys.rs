use std::iter;

/// VCF record samples keys.
#[derive(Debug, Eq, PartialEq)]
pub struct Keys<'a>(&'a str);

impl<'a> Keys<'a> {
    pub(super) fn new(src: &'a str) -> Self {
        Self(src)
    }

    /// Returns an iterator over keys.
    pub fn iter(&self) -> impl Iterator<Item = &'a str> + 'a {
        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(parse_key(&mut src))
            }
        })
    }
}

fn parse_key<'a>(src: &mut &'a str) -> &'a str {
    const DELIMITER: u8 = b':';

    match src.as_bytes().iter().position(|&b| b == DELIMITER) {
        Some(i) => {
            let (buf, rest) = src.split_at(i);
            *src = &rest[1..];
            buf
        }
        None => {
            let (buf, rest) = src.split_at(src.len());
            *src = rest;
            buf
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        let genotypes = Keys::new("");
        assert!(genotypes.iter().next().is_none());

        let keys = Keys::new("GT:GQ");
        let actual: Vec<_> = keys.iter().collect();
        let expected = ["GT", "GQ"];
        assert_eq!(actual, expected);
    }
}
