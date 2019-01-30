use std::str::Split;

const DELIMITER: char  = ';';

pub struct Attributes<'a> {
    split: Split<'a, char>,
}

impl<'a> Attributes<'a> {
    pub fn new(inner: &str) -> Attributes {
        Attributes {
            split: inner.split(DELIMITER),
        }
    }
}

impl<'a> Iterator for Attributes<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        self.split
            .next()
            .map(str::trim_left)
            .and_then(|p| {
                let mut pieces = p.splitn(2, ' ');

                if let Some(key) = pieces.next() {
                    if let Some(value) = pieces.next() {
                        return Some((key, trim_quotes(value)))
                    }
                }

                None
            })
    }
}

fn trim_quotes(s: &str) -> &str {
    s.trim_matches('"')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() {
        let data = r#"gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";"#;
        let mut reader = Attributes::new(data);

        assert_eq!(reader.next(), Some(("gene_id", "ENSG00000223972.5")));
        assert_eq!(reader.next(), Some(("gene_type", "transcribed_unprocessed_pseudogene")));
        assert_eq!(reader.next(), Some(("gene_name", "DDX11L1")));
        assert_eq!(reader.next(), Some(("level", "2")));
        assert_eq!(reader.next(), Some(("havana_gene", "OTTHUMG00000000961.2")));

        assert_eq!(reader.next(), None);
    }

    #[test]
    fn test_trim_quotes() {
        assert_eq!(trim_quotes("DDX11L1"), "DDX11L1");
        assert_eq!(trim_quotes(r#""DDX11L1""#), "DDX11L1");
        assert_eq!(trim_quotes(""), "");
    }
}
