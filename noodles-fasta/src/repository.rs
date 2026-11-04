//! Sequence repository and adapters.

mod adapter;
pub mod adapters;

pub use self::adapter::Adapter;

use std::{
    collections::HashMap,
    fmt, io,
    sync::{Arc, Mutex, RwLock},
};

use noodles_core::region::Interval;

use super::record::Sequence;

// The adapter and the cache have separate locks, so that reading from the adapter (e.g., a range
// read for one slice) does not block cache hits. To avoid deadlocks, the cache lock is never held
// while taking the adapter lock.
struct AdapterCache {
    adapter: Mutex<Box<dyn Adapter>>,
    cache: RwLock<HashMap<Vec<u8>, Arc<Sequence>>>,
}

impl AdapterCache {
    fn cached(&self, name: &[u8]) -> Option<Arc<Sequence>> {
        self.cache.read().unwrap().get(name).cloned()
    }
}

/// A caching sequence repository.
pub struct Repository(Arc<AdapterCache>);

impl Repository {
    /// Creates a sequence repository.
    pub fn new<A>(adapter: A) -> Self
    where
        A: Adapter + 'static,
    {
        Self(Arc::new(AdapterCache {
            adapter: Mutex::new(Box::new(adapter)),
            cache: RwLock::new(HashMap::new()),
        }))
    }

    /// Returns the sequence of the given name.
    pub fn get(&self, name: &[u8]) -> Option<io::Result<Arc<Sequence>>> {
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(sequence));
        }

        let mut adapter = self.0.adapter.lock().unwrap();

        // A concurrent miss may have populated the cache while this call waited.
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(sequence));
        }

        let record = match adapter.get(name)? {
            Ok(record) => record,
            Err(e) => return Some(Err(e)),
        };

        let mut cache = self.0.cache.write().unwrap();

        let s = cache
            .entry(name.into())
            .or_insert_with(|| Arc::new(record.into_sequence()));

        Some(Ok(Arc::clone(s)))
    }

    /// Returns the subsequence in the given interval of the sequence with the given name.
    ///
    /// If the whole sequence is not cached, this reads only the interval, when the adapter
    /// supports it (see [`Adapter::get_range`]), and does not cache it. The subsequence ends
    /// early if the interval extends past the end of the sequence.
    ///
    /// This returns `None` if the sequence is not cached and the adapter does not support
    /// ranges or does not have the sequence. Use [`Self::get`] in that case.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_core::Position;
    /// use noodles_fasta::{self as fasta, record::{Definition, Sequence}, Repository};
    ///
    /// let repository = Repository::new(vec![fasta::Record::new(
    ///     Definition::new("sq0", None),
    ///     Sequence::from(b"ACGT".to_vec()),
    /// )]);
    ///
    /// // The `Vec<Record>` adapter does not read ranges.
    /// let start = Position::try_from(2)?;
    /// let end = Position::try_from(3)?;
    /// assert!(repository.get_range(b"sq0", (start..=end).into()).is_none());
    ///
    /// // Once the whole sequence is cached, its subsequences are, too.
    /// repository.get(b"sq0").transpose()?;
    /// let sequence = repository.get_range(b"sq0", (start..=end).into()).transpose()?;
    /// assert_eq!(sequence.as_deref().map(|s| s.as_ref()), Some(&b"CG"[..]));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get_range(&self, name: &[u8], interval: Interval) -> Option<io::Result<Arc<Sequence>>> {
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(Arc::new(slice_to_end(&sequence, interval))));
        }

        let mut adapter = self.0.adapter.lock().unwrap();

        // The whole sequence may have been loaded while this call waited.
        if let Some(sequence) = self.0.cached(name) {
            return Some(Ok(Arc::new(slice_to_end(&sequence, interval))));
        }

        // Ranges are not cached: a sorted file's slices read mostly disjoint spans, so a later
        // slice would rarely find its span already read.
        Some(adapter.get_range(name, interval)?.map(Arc::new))
    }

    /// Returns the number of cached sequences.
    pub fn len(&self) -> usize {
        self.0.cache.read().unwrap().len()
    }

    /// Returns whether any sequences are cached.
    pub fn is_empty(&self) -> bool {
        self.0.cache.read().unwrap().is_empty()
    }

    /// Clears the sequence cache.
    pub fn clear(&self) {
        self.0.cache.write().unwrap().clear();
    }
}

/// Returns the subsequence in the interval, ending early at the end of the sequence.
fn slice_to_end(sequence: &Sequence, interval: Interval) -> Sequence {
    let src = sequence.as_ref();

    let start = interval.start().map(|p| usize::from(p) - 1).unwrap_or(0);
    let end = interval.end().map(usize::from).unwrap_or(src.len());

    let start = start.min(src.len());
    let end = end.clamp(start, src.len());

    Sequence::from(src[start..end].to_vec())
}

impl Clone for Repository {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl fmt::Debug for Repository {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Repository")
            .field("cache", &*self.0.cache.read().unwrap())
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
    use std::{
        sync::{Barrier, mpsc},
        time::Duration,
    };

    use super::*;
    use crate::{
        Record,
        record::{Definition, Sequence},
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
            Some(Arc::new(sq0.sequence().clone()))
        );
        assert_eq!(repository.get(b"sq1").transpose()?, None);

        Ok(())
    }

    // Loads every name as "ACGT", counting loads. Range reads wait at the barrier when they start
    // and before they load.
    struct TestAdapter(Arc<Mutex<usize>>, Arc<Barrier>);

    impl Adapter for TestAdapter {
        fn get(&mut self, name: &[u8]) -> Option<io::Result<Record>> {
            *self.0.lock().unwrap() += 1;
            std::thread::sleep(Duration::from_millis(10));
            let sequence = Sequence::from(b"ACGT".to_vec());
            Some(Ok(Record::new(Definition::new(name, None), sequence)))
        }

        fn get_range(&mut self, name: &[u8], _: Interval) -> Option<io::Result<Sequence>> {
            self.1.wait();
            self.1.wait();
            Some(self.get(name)?.map(Record::into_sequence))
        }
    }

    #[test]
    fn test_concurrent_cache_misses_load_once() {
        let adapter = TestAdapter(Arc::default(), Arc::new(Barrier::new(2)));
        let loads = Arc::clone(&adapter.0);
        let repository = Repository::new(adapter);
        std::thread::scope(|scope| {
            for _ in 0..4 {
                scope.spawn(|| repository.get(b"sq0").unwrap().unwrap());
            }
        });
        assert_eq!(*loads.lock().unwrap(), 1);
    }

    #[test]
    fn test_cache_hits_do_not_wait_for_adapter_reads() {
        let barrier = Arc::new(Barrier::new(2));
        let repository = &Repository::new(TestAdapter(Arc::default(), Arc::clone(&barrier)));
        repository.get(b"sq0").unwrap().unwrap();
        std::thread::scope(|scope| {
            let range_read = scope.spawn(move || repository.get_range(b"sq1", (..).into()));
            barrier.wait();
            let (tx, rx) = mpsc::channel();
            scope.spawn(move || tx.send(repository.get(b"sq0").is_some()).unwrap());
            let result = rx.recv_timeout(Duration::from_secs(10));
            barrier.wait();
            assert_eq!(result, Ok(true), "a cache hit waited for an adapter read");
            assert!(range_read.join().unwrap().is_some());
        });
    }

    #[test]
    fn test_get_range_matches_whole_sequence() -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Cursor;

        use noodles_core::Region;

        use crate::{fai, io::IndexedReader};

        let data = Cursor::new(b">sq0\nACGT\nACGT\nAC\n>sq1\nNN\n");
        let mut index = fai::io::Reader::new(&b"sq0\t10\t5\t4\t5\nsq1\t2\t23\t2\t3\n"[..]);
        let reader = IndexedReader::new(data, index.read_index()?);
        let ranges = Repository::new(adapters::IndexedReader::new(reader));
        let sequence = Sequence::from(b"ACGTACGTAC".to_vec());
        let whole = Repository::new(vec![Record::new(Definition::new("sq0", None), sequence)]);
        whole.get(b"sq0").transpose()?;

        // Regions and their sequences. Ranges end early at the end of the sequence.
        for case in "sq0:3-6=GTAC sq0:8-12=TAC sq0:6=CGTAC sq0:11= sq0:13-14=".split(' ') {
            let (region, expected) = case.split_once('=').unwrap();
            let region: Region = region.parse()?;
            for r in [&ranges, &whole] {
                let actual = r.get_range(region.name(), region.interval()).transpose()?;
                assert_eq!(
                    actual.as_deref().map(|s| s.as_ref()),
                    Some(expected.as_bytes())
                );
            }
        }

        assert!(ranges.is_empty(), "ranges are not cached");
        assert!(ranges.get_range(b"sqx", (..).into()).unwrap().is_err());
        assert!(whole.get_range(b"sqx", (..).into()).is_none());

        Ok(())
    }
}
