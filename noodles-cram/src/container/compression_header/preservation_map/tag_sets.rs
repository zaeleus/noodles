//! CRAM container preservation map tag sets.

mod key;

pub use self::key::Key;

pub(crate) type TagSets = Vec<Vec<Key>>;
