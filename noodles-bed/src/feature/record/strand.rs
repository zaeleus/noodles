/// A BED record feature strand.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    /// Forward (sense or coding) strand.
    Forward,
    /// Reverse (antisense or complementary) strand.
    Reverse,
}
