use crate::lazy;
use tokio::io::{self, AsyncRead, AsyncReadExt};

pub(super) async fn read_lazy_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    record: &mut lazy::Record,
) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    use crate::reader::record::read_site;

    let l_shared = match reader.read_u32_le().await {
        Ok(n) => usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(0),
        Err(e) => return Err(e),
    };

    let l_indiv = reader.read_u32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    buf.resize(l_shared, Default::default());
    reader.read_exact(buf).await?;
    let mut buf_reader = &buf[..];
    let (n_fmt, n_sample) = read_site(&mut buf_reader, record)?;

    let genotypes = record.genotypes_mut().as_mut();
    genotypes.resize(l_indiv, Default::default());
    reader.read_exact(genotypes).await?;
    record.genotypes_mut().set_format_count(n_fmt);
    record.genotypes_mut().set_sample_count(n_sample);

    Ok(l_shared + l_indiv)
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        header::{format, info},
        record::{
            genotypes::{sample::Value as GenotypeFieldValue, Keys},
            info::field::Value as InfoFieldValue,
            Filters as VcfFilters, Genotypes as VcfGenotypes, Ids, Position, QualityScore,
        },
    };

    use super::*;
    use crate::header::StringMaps;

    #[tokio::test]
    async fn test_read_lazy_record() -> Result<(), Box<dyn std::error::Error>> {
        use crate::reader::record::tests::{DATA, RAW_HEADER};

        let header = RAW_HEADER.parse()?;
        let string_maps: StringMaps = RAW_HEADER.parse()?;

        let mut reader = &DATA[..];
        let mut buf = Vec::new();
        let mut record = lazy::Record::default();
        read_lazy_record(&mut reader, &mut buf, &mut record).await?;

        assert_eq!(record.chromosome_id(), 1);
        assert_eq!(record.position(), Position::from(101));
        assert_eq!(record.rlen(), 1);
        assert_eq!(
            record.quality_score(),
            QualityScore::try_from(30.1).map(Some)?
        );
        assert_eq!(record.ids(), &"rs123".parse::<Ids>()?);
        assert_eq!(record.reference_bases(), &"A".parse()?);
        assert_eq!(record.alternate_bases(), &"C".parse()?);

        assert_eq!(
            record
                .filters()
                .try_into_vcf_record_filters(string_maps.strings())?,
            Some(VcfFilters::Pass),
        );

        // info

        let actual = record
            .info()
            .try_into_vcf_record_info(&header, string_maps.strings())?;

        let expected = [
            ("HM3".parse()?, Some(InfoFieldValue::Flag)),
            (info::key::ALLELE_COUNT, Some(InfoFieldValue::Integer(3))),
            (
                info::key::TOTAL_ALLELE_COUNT,
                Some(InfoFieldValue::Integer(6)),
            ),
            (
                info::key::ANCESTRAL_ALLELE,
                Some(InfoFieldValue::String(String::from("C"))),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);

        // genotypes

        let actual = record
            .genotypes()
            .try_into_vcf_record_genotypes(&header, string_maps.strings())?;

        let expected = VcfGenotypes::new(
            Keys::try_from(vec![
                format::key::GENOTYPE,
                format::key::CONDITIONAL_GENOTYPE_QUALITY,
                format::key::READ_DEPTH,
                format::key::READ_DEPTHS,
                format::key::ROUNDED_GENOTYPE_LIKELIHOODS,
            ])?,
            vec![
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/0"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(32)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(32), Some(0)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(0),
                        Some(10),
                        Some(100),
                    ])),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("0/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(48)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(32), Some(16)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(10),
                        Some(0),
                        Some(100),
                    ])),
                ],
                vec![
                    Some(GenotypeFieldValue::String(String::from("1/1"))),
                    Some(GenotypeFieldValue::Integer(10)),
                    Some(GenotypeFieldValue::Integer(64)),
                    Some(GenotypeFieldValue::IntegerArray(vec![Some(0), Some(64)])),
                    Some(GenotypeFieldValue::IntegerArray(vec![
                        Some(100),
                        Some(10),
                        Some(0),
                    ])),
                ],
            ],
        );

        assert_eq!(actual, expected);

        Ok(())
    }
}
