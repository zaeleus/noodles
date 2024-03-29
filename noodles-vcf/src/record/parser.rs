use super::Record;
use crate::{reader::record::ParseError, Header};

pub fn parse(s: &str, header: &Header) -> Result<Record, ParseError> {
    use crate::reader::parse_record;

    let mut record = Record::default();
    parse_record(s, header, &mut record)?;
    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::{alternate_bases::Allele, reference_bases::Base, Chromosome, Filters};

        let header = Header::default();
        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL";
        let record = parse(s, &header)?;

        assert!(matches!(record.chromosome(), Chromosome::Name(name) if name == "chr1"));

        assert_eq!(usize::from(record.position()), 13);

        let ids = "nd0".parse()?;
        assert_eq!(record.ids(), &ids);

        let reference_bases = [Base::A, Base::T, Base::C, Base::G];
        assert_eq!(&record.reference_bases()[..], &reference_bases[..]);

        let alternate_bases = [Allele::Bases(vec![Base::A])];
        assert_eq!(&record.alternate_bases()[..], &alternate_bases[..]);

        assert_eq!(record.quality_score().map(f32::from), Some(5.8));
        assert_eq!(record.filters(), Some(&Filters::Pass));
        assert_eq!(record.info().len(), 1);
        assert!(record.genotypes().is_empty());

        Ok(())
    }

    #[test]
    fn test_parse_with_genotype_info() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::{
            genotypes::{keys::key, sample::Value, Keys},
            Genotypes,
        };

        let header = Header::builder().add_sample_name("sample0").build();
        let s = "chr1\t13\tnd0\tATCG\tA\t5.8\tPASS\tSVTYPE=DEL\tGT:GQ\t0|1:13";
        let record = parse(s, &header)?;

        let keys = Keys::try_from(vec![key::GENOTYPE, key::CONDITIONAL_GENOTYPE_QUALITY])?;
        let values = vec![vec![
            Some(Value::String(String::from("0|1"))),
            Some(Value::Integer(13)),
        ]];

        let actual = record.genotypes();
        let expected = Genotypes::new(keys, values);

        assert_eq!(actual, &expected);

        Ok(())
    }
}
