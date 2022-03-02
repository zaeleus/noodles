//! Sequence repository and adapters.

mod adapter;
pub mod adapters;

pub use self::adapter::Adapter;

use std::{
    collections::HashMap,
    io,
    sync::{Arc, RwLock},
};

use super::record::Sequence;

#[derive(Debug)]
struct AdapterCache<A> {
    adapter: A,
    cache: HashMap<String, Sequence>,
}

/// A caching sequence repository.
#[derive(Debug)]
pub struct Repository<A>(Arc<RwLock<AdapterCache<A>>>);

impl<A> Repository<A>
where
    A: Adapter,
{
    /// Creates a sequence repository.
    pub fn new(adapter: A) -> Self {
        Self(Arc::new(RwLock::new(AdapterCache {
            adapter,
            cache: HashMap::new(),
        })))
    }

    /// Returns the sequence of the given name.
    pub fn get(&self, name: &str) -> Option<io::Result<Sequence>> {
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

impl<A> Clone for Repository<A>
where
    A: Adapter,
{
    fn clone(&self) -> Self {
        Self(self.0.clone())
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
            repository.get("sq0").transpose()?,
            Some(sq0.sequence().clone())
        );
        assert_eq!(repository.get("sq1").transpose()?, None);

        Ok(())
    }
}
