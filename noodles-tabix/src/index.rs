//! Tabix index.

mod indexer;

pub use self::indexer::Indexer;

pub(crate) const DEPTH: u8 = 5;
