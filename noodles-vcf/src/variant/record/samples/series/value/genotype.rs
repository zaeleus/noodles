//! Variant record samples genotype value.

mod phasing;

use std::{fmt::Debug, io};

pub use self::phasing::Phasing;

/// A variant record samples genotype value.
pub trait Genotype: Debug {
    /// Returns an iterator over allele position-phasing pairs.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_>;
}
