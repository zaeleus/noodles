#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Version {
    major: u8,
    minor: u8,
}

impl Version {
    pub fn new(major: u8, minor: u8) -> Self {
        Self { major, minor }
    }

    pub fn major(&self) -> u8 {
        self.major
    }

    pub fn minor(&self) -> u8 {
        self.minor
    }
}

impl Default for Version {
    fn default() -> Self {
        Self::new(3, 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Version::default(), Version::new(3, 0));
    }
}
