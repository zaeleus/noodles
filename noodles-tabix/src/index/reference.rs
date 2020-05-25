pub mod bin;

pub use self::bin::Bin;

#[derive(Debug)]
pub struct Reference {
    bins: Vec<Bin>,
    intervals: Vec<u64>,
}

impl Reference {
    pub fn new(bins: Vec<Bin>, intervals: Vec<u64>) -> Self {
        Self { bins, intervals }
    }
}
