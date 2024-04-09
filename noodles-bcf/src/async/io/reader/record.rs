use std::mem;

use tokio::io::{self, AsyncRead, AsyncReadExt};

use crate::Record;

pub(super) async fn read_record<R>(reader: &mut R, record: &mut Record) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let l_shared = read_site_length(reader).await?;
    let l_indiv = read_samples_length(reader).await?;

    let site_buf = record.fields_mut().site_buf_mut();
    site_buf.resize(l_shared, 0);
    reader.read_exact(site_buf).await?;

    record.fields_mut().index()?;

    let samples_buf = record.fields_mut().samples_buf_mut();
    samples_buf.resize(l_indiv, 0);
    reader.read_exact(samples_buf).await?;

    Ok(l_shared + l_indiv)
}

async fn read_site_length<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    let mut buf = [0; mem::size_of::<u32>()];
    read_exact_or_eof(reader, &mut buf).await?;
    let n = u32::from_le_bytes(buf);
    usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_exact_or_eof<R>(reader: &mut R, mut buf: &mut [u8]) -> io::Result<()>
where
    R: AsyncRead + Unpin,
{
    let mut bytes_read = 0;

    while !buf.is_empty() {
        match reader.read(buf).await {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                bytes_read += n;
            }
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    if bytes_read > 0 && !buf.is_empty() {
        Err(io::Error::new(io::ErrorKind::UnexpectedEof, "early eof"))
    } else {
        Ok(())
    }
}

async fn read_samples_length<R>(reader: &mut R) -> io::Result<usize>
where
    R: AsyncRead + Unpin,
{
    reader
        .read_u32_le()
        .await
        .and_then(|n| usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;
    use noodles_vcf::{
        self as vcf,
        variant::{
            record::{
                info::{self, field::Value as InfoFieldValue},
                samples::{self, series::value::genotype::Phasing, Sample},
                AlternateBases, Filters, Info, Samples,
            },
            record_buf::{
                samples::sample::{value::genotype::Allele, Value as GenotypeFieldValue},
                Samples as VcfGenotypes,
            },
        },
    };

    use super::*;

    #[tokio::test]
    async fn test_read_record() -> Result<(), Box<dyn std::error::Error>> {
        use crate::io::reader::record::tests::{DATA, RAW_HEADER};

        let mut header: vcf::Header = RAW_HEADER.parse()?;
        *header.string_maps_mut() = RAW_HEADER.parse()?;

        let mut reader = &DATA[..];
        let mut record = Record::default();
        read_record(&mut reader, &mut record).await?;

        assert_eq!(record.reference_sequence_id()?, 1);
        assert_eq!(
            record.variant_start().transpose()?,
            Some(Position::try_from(101)?)
        );
        assert_eq!(record.rlen()?, 1);
        assert_eq!(record.quality_score()?, Some(30.1));
        assert_eq!(record.ids().as_ref(), b"rs123");
        assert_eq!(record.reference_bases().as_ref(), b"A");
        assert_eq!(
            record
                .alternate_bases()
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            ["C"]
        );

        assert_eq!(
            record
                .filters()
                .iter(&header)
                .collect::<io::Result<Vec<_>>>()?,
            ["PASS"],
        );

        // info

        let info = record.info();
        let mut iter = info.iter(&header);
        assert!(matches!(
            iter.next(),
            Some(Ok(("HM3", Some(InfoFieldValue::Flag))))
        ));
        assert!(matches!(
            iter.next(),
            Some(Ok((
                info::field::key::ALLELE_COUNT,
                Some(InfoFieldValue::Integer(3))
            )))
        ));
        assert!(matches!(
            iter.next(),
            Some(Ok((
                info::field::key::TOTAL_ALLELE_COUNT,
                Some(InfoFieldValue::Integer(6))
            )))
        ));
        assert!(matches!(
            iter.next(),
            Some(Ok((
                info::field::key::ANCESTRAL_ALLELE,
                Some(InfoFieldValue::String("C"))
            )))
        ));
        assert!(iter.next().is_none());

        // genotypes

        let samples = record.samples()?;

        let keys = samples
            .column_names(&header)
            .map(|result| result.map(String::from))
            .collect::<io::Result<_>>()?;

        let values = samples
            .iter()
            .map(|sample| {
                sample
                    .iter(&header)
                    .map(|result| {
                        result.and_then(|(_, value)| value.map(|v| v.try_into()).transpose())
                    })
                    .collect()
            })
            .collect::<io::Result<_>>()?;

        let actual = VcfGenotypes::new(keys, values);

        let expected = VcfGenotypes::new(
            [
                String::from(samples::keys::key::GENOTYPE),
                String::from(samples::keys::key::CONDITIONAL_GENOTYPE_QUALITY),
                String::from(samples::keys::key::READ_DEPTH),
                String::from(samples::keys::key::READ_DEPTHS),
                String::from(samples::keys::key::ROUNDED_GENOTYPE_LIKELIHOODS),
            ]
            .into_iter()
            .collect(),
            vec![
                vec![
                    Some(GenotypeFieldValue::Genotype(
                        [
                            Allele::new(Some(0), Phasing::Unphased),
                            Allele::new(Some(0), Phasing::Unphased),
                        ]
                        .into_iter()
                        .collect(),
                    )),
                    Some(GenotypeFieldValue::from(10)),
                    Some(GenotypeFieldValue::from(32)),
                    Some(GenotypeFieldValue::from(vec![Some(32), Some(0)])),
                    Some(GenotypeFieldValue::from(vec![Some(0), Some(10), Some(100)])),
                ],
                vec![
                    Some(GenotypeFieldValue::Genotype(
                        [
                            Allele::new(Some(0), Phasing::Unphased),
                            Allele::new(Some(1), Phasing::Unphased),
                        ]
                        .into_iter()
                        .collect(),
                    )),
                    Some(GenotypeFieldValue::from(10)),
                    Some(GenotypeFieldValue::from(48)),
                    Some(GenotypeFieldValue::from(vec![Some(32), Some(16)])),
                    Some(GenotypeFieldValue::from(vec![Some(10), Some(0), Some(100)])),
                ],
                vec![
                    Some(GenotypeFieldValue::Genotype(
                        [
                            Allele::new(Some(1), Phasing::Unphased),
                            Allele::new(Some(1), Phasing::Unphased),
                        ]
                        .into_iter()
                        .collect(),
                    )),
                    Some(GenotypeFieldValue::from(10)),
                    Some(GenotypeFieldValue::from(64)),
                    Some(GenotypeFieldValue::from(vec![Some(0), Some(64)])),
                    Some(GenotypeFieldValue::from(vec![Some(100), Some(10), Some(0)])),
                ],
            ],
        );

        assert_eq!(actual, expected);

        Ok(())
    }

    #[tokio::test]
    async fn test_read_site_length() -> io::Result<()> {
        let data = [0x08, 0x00, 0x00, 0x00];
        let mut reader = &data[..];
        assert_eq!(read_site_length(&mut reader).await?, 8);

        let data = [];
        let mut reader = &data[..];
        assert_eq!(read_site_length(&mut reader).await?, 0);

        let data = [0x08];
        let mut reader = &data[..];
        assert!(matches!(
            read_site_length(&mut reader).await,
            Err(e) if e.kind() == io::ErrorKind::UnexpectedEof
        ));

        Ok(())
    }
}
