/// A VCF record genotype field array value.
#[derive(Clone, Debug, PartialEq)]
pub enum Array {
    /// An array of 32-bit integers.
    Integer(Vec<Option<i32>>),
    /// An array of single-precision floating-points.
    Float(Vec<Option<f32>>),
    /// An array of characters.
    Character(Vec<Option<char>>),
    /// An array of strings.
    String(Vec<Option<String>>),
}
