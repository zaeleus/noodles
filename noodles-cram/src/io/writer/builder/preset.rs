#[derive(Default)]
pub enum Preset {
    #[default]
    Medium,
}

impl Preset {
    pub fn records_per_slice(&self) -> usize {
        match self {
            Self::Medium => 10240,
        }
    }
}
