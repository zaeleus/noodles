/// A variant record samples series genotype value phasing.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Phasing {
    /// Phased.
    Phased,
    /// Unphased.
    Unphased,
}
