pub(super) fn parse_reference_sequence_name(s: &str, reference_sequence_name: &mut String) {
    reference_sequence_name.clear();
    reference_sequence_name.push_str(s);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_reference_sequence_name() {
        let mut buf = String::new();

        buf.clear();
        let src = "sq0";
        parse_reference_sequence_name(src, &mut buf);
        assert_eq!(buf, "sq0");

        buf.clear();
        let src = "<sq0>";
        parse_reference_sequence_name(src, &mut buf);
        assert_eq!(buf, "<sq0>");
    }
}
