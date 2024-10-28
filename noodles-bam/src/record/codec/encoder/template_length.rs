use super::num::write_i32_le;

pub(super) fn write_template_length(dst: &mut Vec<u8>, template_length: i32) {
    write_i32_le(dst, template_length);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_template_length() {
        let mut buf = Vec::new();
        write_template_length(&mut buf, 8);
        assert_eq!(buf, [0x08, 0x00, 0x00, 0x00]);
    }
}
