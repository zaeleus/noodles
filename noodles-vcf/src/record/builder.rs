use std::{error, fmt};

use super::{
    AlternateBases, Chromosome, FilterStatus, Format, Genotype, Id, Info, QualityScore, Record,
    ReferenceBases,
};

#[derive(Debug, Default)]
pub struct Builder {
    chromosome: Option<Chromosome>,
    position: Option<i32>,
    id: Id,
    reference_bases: ReferenceBases,
    alternate_bases: AlternateBases,
    quality_score: QualityScore,
    filter_status: FilterStatus,
    info: Info,
    format: Option<Format>,
    genotypes: Vec<Genotype>,
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
        self.position = Some(position);
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

    pub fn set_quality_score(mut self, quality_score: QualityScore) -> Self {
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

    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    pub fn set_genotypes(mut self, genotypes: Vec<Genotype>) -> Self {
        self.genotypes = genotypes;
        self
    }

    pub fn build(self) -> Result<Record, BuildError> {
        if self.reference_bases.is_empty() {
            return Err(BuildError);
        }

        Ok(Record {
            chromosome: self.chromosome.ok_or_else(|| BuildError)?,
            position: self.position.ok_or_else(|| BuildError)?,
            id: self.id,
            reference_bases: self.reference_bases,
            alternate_bases: self.alternate_bases,
            quality_score: self.quality_score,
            filter_status: self.filter_status,
            info: self.info,
            format: self.format,
            genotypes: self.genotypes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() -> Result<(), Box<dyn std::error::Error>> {
        let chromosome: Chromosome = "sq0".parse()?;
        let reference_bases: ReferenceBases = "A".parse()?;

        let record = Builder::new()
            .set_chromosome(chromosome.clone())
            .set_position(5)
            .set_reference_bases(reference_bases.clone())
            .build()?;

        assert_eq!(record.chromosome(), &chromosome);
        assert_eq!(record.position(), 5);
        assert!(record.id().is_none());
        assert_eq!(record.reference_bases(), &reference_bases);
        assert!(record.alternate_bases().is_empty());
        assert!(record.quality_score().is_none());
        assert_eq!(record.filter_status(), &FilterStatus::Missing);
        assert!(record.info().is_empty());
        assert!(record.format().is_none());
        assert!(record.genotypes().is_empty());

        Ok(())
    }
}
