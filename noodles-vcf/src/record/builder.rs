use std::{convert::TryFrom, error, fmt};

use super::{
    reference_bases::Base, AlternateBases, Chromosome, FilterStatus, Format, Genotype, Ids, Info,
    QualityScore, Record, ReferenceBases,
};

#[derive(Debug, Default, PartialEq)]
pub struct Builder {
    chromosome: Option<Chromosome>,
    position: Option<i32>,
    ids: Ids,
    reference_bases: Vec<Base>,
    alternate_bases: AlternateBases,
    quality_score: QualityScore,
    filter_status: FilterStatus,
    info: Info,
    format: Option<Format>,
    genotypes: Vec<Genotype>,
}

/// An error returned when a VCF record fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The chromosome is missing.
    MissingChromosome,
    /// The position is missing.
    MissingPosition,
    /// The reference bases are missing.
    MissingReferenceBases,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingChromosome => f.write_str("missing chromosome"),
            Self::MissingPosition => f.write_str("missing position"),
            Self::MissingReferenceBases => f.write_str("missing reference bases"),
        }
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

    pub fn set_ids(mut self, ids: Ids) -> Self {
        self.ids = ids;
        self
    }

    pub fn set_reference_bases(mut self, reference_bases: ReferenceBases) -> Self {
        self.reference_bases.clear();
        self.reference_bases.extend(reference_bases.iter());
        self
    }

    pub fn add_reference_base(mut self, reference_base: Base) -> Self {
        self.reference_bases.push(reference_base);
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

    pub fn add_genotype(mut self, genotype: Genotype) -> Self {
        self.genotypes.push(genotype);
        self
    }

    pub fn build(self) -> Result<Record, BuildError> {
        Ok(Record {
            chromosome: self
                .chromosome
                .ok_or_else(|| BuildError::MissingChromosome)?,
            position: self.position.ok_or_else(|| BuildError::MissingPosition)?,
            ids: self.ids,
            reference_bases: ReferenceBases::try_from(self.reference_bases)
                .map_err(|_| BuildError::MissingReferenceBases)?,
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
        assert!(record.ids().is_empty());
        assert_eq!(record.reference_bases(), &reference_bases);
        assert!(record.alternate_bases().is_empty());
        assert!(record.quality_score().is_none());
        assert_eq!(record.filter_status(), &FilterStatus::Missing);
        assert!(record.info().is_empty());
        assert!(record.format().is_none());
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_build() -> Result<(), Box<dyn std::error::Error>> {
        let result = Builder::new()
            .set_position(1)
            .set_reference_bases("A".parse()?)
            .build();
        assert_eq!(result, Err(BuildError::MissingChromosome));

        let result = Builder::new()
            .set_chromosome("sq0".parse()?)
            .set_reference_bases("A".parse()?)
            .build();
        assert_eq!(result, Err(BuildError::MissingPosition));

        let result = Builder::new()
            .set_chromosome("sq0".parse()?)
            .set_position(1)
            .build();
        assert_eq!(result, Err(BuildError::MissingReferenceBases));

        Ok(())
    }
}
