pub(super) fn parse_chromosome(s: &str, chromosome: &mut String) {
    chromosome.clear();
    chromosome.push_str(s);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_chromosome() {
        let mut buf = String::new();

        buf.clear();
        let src = "sq0";
        parse_chromosome(src, &mut buf);
        assert_eq!(buf, "sq0");

        buf.clear();
        let src = "<sq0>";
        parse_chromosome(src, &mut buf);
        assert_eq!(buf, "<sq0>");
    }
}
