use super::Value;

#[derive(Debug)]
pub struct Field {
    tag: String,
    value: Value,
}

impl Field {
    pub fn new<T>(tag: T, value: Value) -> Field where T: Into<String> {
        Field { tag: tag.into(), value }
    }

    pub fn tag(&self) -> &str {
        &self.tag
    }

    pub fn value(&self) -> &Value {
        &self.value
    }
}
