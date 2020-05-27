pub mod reference;

pub use self::reference::Reference;

#[derive(Debug)]
pub struct Index {
    references: Vec<Reference>,
    n_no_coor: Option<u64>,
}

impl Index {
    pub fn new(references: Vec<Reference>, n_no_coor: Option<u64>) -> Self {
        Self {
            references,
            n_no_coor,
        }
    }

    pub fn references(&self) -> &[Reference] {
        &self.references
    }

    pub fn unplaced_unmapped_read_count(&self) -> Option<u64> {
        self.n_no_coor
    }
}
