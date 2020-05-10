use std::{error, fmt};

use super::{AlternateBases, Chromosome, FilterStatus, Id, Info, Record, ReferenceBases};

#[derive(Debug, Default)]
pub struct Builder {
    chromosome: Option<Chromosome>,
    position: i32,
    id: Id,
    reference_bases: ReferenceBases,
    alternate_bases: AlternateBases,
    quality_score: f32,
    filter_status: FilterStatus,
    info: Info,
}

#[derive(Debug)]
pub struct BuildError;

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid builder state")
    }
}
impl Builder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_chromosome(mut self, chromosome: Chromosome) -> Self {
        self.chromosome = Some(chromosome);
        self
    }

    pub fn set_position(mut self, position: i32) -> Self {
        self.position = position;
        self
    }

    pub fn set_id(mut self, id: Id) -> Self {
        self.id = id;
        self
    }

    pub fn set_reference_bases(mut self, reference_bases: ReferenceBases) -> Self {
        self.reference_bases = reference_bases;
        self
    }

    pub fn set_alternate_bases(mut self, alternate_bases: AlternateBases) -> Self {
        self.alternate_bases = alternate_bases;
        self
    }

    pub fn set_quality_score(mut self, quality_score: f32) -> Self {
        self.quality_score = quality_score;
        self
    }

    pub fn set_filter_status(mut self, filter_status: FilterStatus) -> Self {
        self.filter_status = filter_status;
        self
    }

    pub fn set_info(mut self, info: Info) -> Self {
        self.info = info;
        self
    }

    pub fn build(self) -> Result<Record, BuildError> {
        Ok(Record {
            chromosome: self.chromosome.ok_or_else(|| BuildError)?,
            position: self.position,
            id: self.id,
            reference_bases: self.reference_bases,
            alternate_bases: self.alternate_bases,
            quality_score: self.quality_score,
            filter_status: self.filter_status,
            info: self.info,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let chromosome: Chromosome = "sq0".parse()?;
        let record = Builder::new().set_chromosome(chromosome.clone()).build()?;

        assert_eq!(record.chromosome(), &chromosome);
        assert_eq!(record.position(), 0);
        assert!(record.id().is_none());
        assert!(record.reference_bases().is_empty());
        assert!(record.alternate_bases().is_empty());
        assert_eq!(record.quality_score(), 0.0);
        assert_eq!(record.filter_status(), &FilterStatus::Missing);
        assert!(record.info().is_empty());

        Ok(())
    }
}
