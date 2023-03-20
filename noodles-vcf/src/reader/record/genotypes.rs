mod keys;
mod values;

use std::io;

use self::{keys::parse_keys, values::parse_values};
use super::next_field;
use crate::{record::Genotypes, Header};

pub(super) fn parse_genotypes(
    header: &Header,
    mut s: &str,
    genotypes: &mut Genotypes,
) -> io::Result<()> {
    genotypes.keys.clear();

    for values in &mut genotypes.values {
        values.clear();
    }

    if s.is_empty() {
        return Ok(());
    }

    let field = next_field(&mut s);
    parse_keys(header, field, &mut genotypes.keys)?;

    genotypes
        .values
        .resize(header.sample_names().len(), Vec::new());

    for values in &mut genotypes.values {
        let field = next_field(&mut s);
        parse_values(header, &genotypes.keys, field, values)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genotypes() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header::format::key,
            record::genotypes::{sample::Value, Keys},
        };

        let mut genotypes = Genotypes::default();

        let header = Header::default();
        parse_genotypes(&header, "", &mut genotypes)?;
        assert!(genotypes.is_empty());

        let header = Header::builder().add_sample_name("sample0").build();
        parse_genotypes(&header, "GT\t0|0", &mut genotypes)?;
        let expected = Genotypes::new(
            Keys::try_from(vec![key::GENOTYPE])?,
            vec![vec![Some(Value::String(String::from("0|0")))]],
        );
        assert_eq!(genotypes, expected);

        let header = Header::builder()
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();
        parse_genotypes(&header, "GQ\t8\t13", &mut genotypes)?;
        let expected = Genotypes::new(
            Keys::try_from(vec![key::CONDITIONAL_GENOTYPE_QUALITY])?,
            vec![
                vec![Some(Value::Integer(8))],
                vec![Some(Value::Integer(13))],
            ],
        );
        assert_eq!(genotypes, expected);

        let header = Header::builder().add_sample_name("sample0").build();
        assert!(matches!(
            parse_genotypes(&header, "GT:GQ", &mut genotypes),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_genotypes(&header, "\t0|0", &mut genotypes),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_genotypes(&header, "GQ\tndls", &mut genotypes),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
