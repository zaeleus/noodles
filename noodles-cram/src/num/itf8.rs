use super::Itf8;

pub fn size_of(value: Itf8) -> usize {
    if value >> (8 - 1) == 0 {
        1
    } else if value >> (16 - 2) == 0 {
        2
    } else if value >> (24 - 3) == 0 {
        3
    } else if value >> (32 - 4) == 0 {
        4
    } else {
        5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_size_of() {
        assert_eq!(size_of(0), 1);
        assert_eq!(size_of(1877), 2);
        assert_eq!(size_of(480665), 3);
        assert_eq!(size_of(123050342), 4);
        assert_eq!(size_of(1968805474), 5);
        assert_eq!(size_of(-1), 5);
    }
}
