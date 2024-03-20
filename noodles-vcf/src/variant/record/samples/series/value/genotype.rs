use std::io;

/// A variant record samples genotype value.
pub trait Genotype {
    /// Returns an iterator over allele position-phasing pairs.
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, u8)>> + '_>;
}
