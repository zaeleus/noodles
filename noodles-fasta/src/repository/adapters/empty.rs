use std::io;

use crate::{repository::Adapter, Record};

/// An empty adapter.
///
/// This adapter always returns `None`.
#[derive(Default)]
pub struct Empty;

impl Empty {
    /// Creates an empty adapter.
    pub fn new() -> Self {
        Self
    }
}

impl Adapter for Empty {
    fn get(&mut self, _: &[u8]) -> Option<io::Result<Record>> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get() {
        let mut adapter = Empty::new();
        assert!(adapter.get(b"sq0").is_none());
    }
}
