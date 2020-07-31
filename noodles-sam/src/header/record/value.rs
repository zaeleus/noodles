/// A SAM header record value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Value {
    /// A string.
    String(String),
    /// A list of key-value pairs.
    Map(Vec<(String, String)>),
}
