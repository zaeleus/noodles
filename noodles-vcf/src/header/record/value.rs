/// A VCF header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A structure.
    Struct(Vec<(String, String)>),
}
