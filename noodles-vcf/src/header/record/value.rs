#[derive(Debug, Eq, PartialEq)]
pub enum Value {
    String(String),
    Struct(Vec<(String, String)>),
}
