//! Sequence repository and adapters.

mod adapter;
pub mod adapters;

pub use self::adapter::Adapter;

use std::{
    collections::HashMap,
    fmt, io,
    sync::{Arc, RwLock},
};

use super::record::Sequence;

struct AdapterCache {
    adapter: Box<dyn Adapter>,
    cache: HashMap<Vec<u8>, Sequence>,
}

/// A caching sequence repository.
pub struct Repository(Arc<RwLock<AdapterCache>>);

impl Repository {
    /// Creates a sequence repository.
    pub fn new<A>(adapter: A) -> Self
    where
        A: Adapter + 'static,
    {
        Self(Arc::new(RwLock::new(AdapterCache {
            adapter: Box::new(adapter),
            cache: HashMap::new(),
        })))
    }

    /// Returns the sequence of the given name.
    pub fn get(&self, name: &[u8]) -> Option<io::Result<Sequence>> {
        {
            let lock = self.0.read().unwrap();

            if let Some(sequence) = lock.cache.get(name) {
                return Some(Ok(sequence.clone()));
            }
        }

        let mut lock = self.0.write().unwrap();

        let record = match lock.adapter.get(name)? {
            Ok(record) => record,
            Err(e) => return Some(Err(e)),
        };

        lock.cache
            .entry(name.into())
            .or_insert_with(|| record.sequence().clone());

        Some(Ok(record.sequence().clone()))
    }

    /// Returns the number of cached sequences.
    pub fn len(&self) -> usize {
        self.0.read().unwrap().cache.len()
    }

    /// Returns whether any sequences are cached.
    pub fn is_empty(&self) -> bool {
        self.0.read().unwrap().cache.is_empty()
    }

    /// Clears the sequence cache.
    pub fn clear(&self) {
        self.0.write().unwrap().cache.clear();
    }
}

impl Clone for Repository {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl fmt::Debug for Repository {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Repository")
            .field("cache", &self.0.read().unwrap().cache)
            .finish()
    }
}

impl Default for Repository {
    fn default() -> Self {
        Self::new(adapters::Empty::new())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        record::{Definition, Sequence},
        Record,
    };

    #[test]
    fn test_get() -> io::Result<()> {
        let sq0 = Record::new(
            Definition::new("sq0", None),
            Sequence::from(b"ACGT".to_vec()),
        );
        let repository = Repository::new(vec![sq0.clone()]);

        assert_eq!(
            repository.get(b"sq0").transpose()?,
            Some(sq0.sequence().clone())
        );
        assert_eq!(repository.get(b"sq1").transpose()?, None);

        Ok(())
    }
}
