/// A description on how skipped sequence bases should be interpreted.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Status {
    /// Low probability of modification (`.`).
    #[default]
    Implicit,
    /// No information (`?`).
    Explicit,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(Status::default(), Status::Implicit);
    }
}
