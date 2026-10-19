use std::io::{self, Write};

use crate::{
    Header,
    header::FileFormat,
    variant::record::samples::series::value::{Genotype, genotype::Phasing},
};

const VCF_4_4: FileFormat = FileFormat::new(4, 4);

pub(super) fn write_genotype<W>(
    writer: &mut W,
    header: &Header,
    genotype: &dyn Genotype,
) -> io::Result<()>
where
    W: Write,
{
    let alleles: Vec<_> = genotype.iter().collect::<io::Result<_>>()?;
    write_alleles(writer, header.file_format(), &alleles)
}

fn write_alleles<W>(
    writer: &mut W,
    file_format: FileFormat,
    alleles: &[(Option<usize>, Phasing)],
) -> io::Result<()>
where
    W: Write,
{
    let (&(first_position, first_phasing), rest) = alleles
        .split_first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "empty input"))?;

    let is_implicitly_unphased = rest
        .iter()
        .any(|&(_, phasing)| phasing == Phasing::Unphased);

    let implicit_phasing = if is_implicitly_unphased {
        Phasing::Unphased
    } else {
        Phasing::Phased
    };

    if file_format < VCF_4_4 {
        if first_phasing != implicit_phasing {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "invalid phasing",
            ));
        }
    } else if first_phasing != implicit_phasing {
        write_phasing(writer, first_phasing)?;
    }

    write_position(writer, first_position)?;

    for &(position, phasing) in rest {
        write_phasing(writer, phasing)?;
        write_position(writer, position)?;
    }

    Ok(())
}

fn write_phasing<W>(writer: &mut W, phasing: Phasing) -> io::Result<()>
where
    W: Write,
{
    const PHASED: u8 = b'|';
    const UNPHASED: u8 = b'/';

    match phasing {
        Phasing::Phased => writer.write_all(&[PHASED]),
        Phasing::Unphased => writer.write_all(&[UNPHASED]),
    }
}

fn write_position<W>(writer: &mut W, position: Option<usize>) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = b'.';

    if let Some(n) = position {
        write!(writer, "{n}")
    } else {
        writer.write_all(&[MISSING])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record_buf::samples::sample::value::{
        Genotype as GenotypeBuf, genotype::Allele,
    };

    #[test]
    fn test_vcf_4_0_write_genotype() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .set_file_format(FileFormat::new(4, 3))
            .build();

        let mut buf = Vec::new();

        buf.clear();
        let genotype = &[Allele::new(Some(0), Phasing::Phased)]
            .into_iter()
            .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0");

        buf.clear();
        let genotype = &[Allele::new(Some(0), Phasing::Unphased)]
            .into_iter()
            .collect();
        assert!(matches!(
            write_genotype(&mut buf, &header, &genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0/1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0|1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        assert!(matches!(
            write_genotype(&mut buf, &header, &genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        assert!(matches!(
            write_genotype(&mut buf, &header, &genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0/1|2");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        assert!(matches!(
            write_genotype(&mut buf, &header, &genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let genotype = &[
            Allele::new(None, Phasing::Unphased),
            Allele::new(None, Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"./.");

        buf.clear();
        let genotype = GenotypeBuf::default();
        assert!(matches!(
            write_genotype(&mut buf, &header, &&genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_vcf_4_4_write_genotype() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .set_file_format(FileFormat::new(4, 4))
            .build();

        let mut buf = Vec::new();

        buf.clear();
        let genotype = &[Allele::new(Some(0), Phasing::Phased)]
            .into_iter()
            .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0");

        buf.clear();
        let genotype = &[Allele::new(Some(0), Phasing::Unphased)]
            .into_iter()
            .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"/0");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0/1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0|1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"/0|1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"|0/1");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"0/1|2");

        buf.clear();
        let genotype = &[
            Allele::new(Some(0), Phasing::Phased),
            Allele::new(Some(1), Phasing::Unphased),
            Allele::new(Some(2), Phasing::Phased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"|0/1|2");

        buf.clear();
        let genotype = &[
            Allele::new(None, Phasing::Unphased),
            Allele::new(None, Phasing::Unphased),
        ]
        .into_iter()
        .collect();
        write_genotype(&mut buf, &header, &genotype)?;
        assert_eq!(buf, b"./.");

        buf.clear();
        let genotype = GenotypeBuf::default();
        assert!(matches!(
            write_genotype(&mut buf, &header, &&genotype),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_phasing() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_phasing(&mut buf, Phasing::Phased)?;
        assert_eq!(buf, b"|");

        buf.clear();
        write_phasing(&mut buf, Phasing::Unphased)?;
        assert_eq!(buf, b"/");

        Ok(())
    }

    #[test]
    fn test_write_position() -> io::Result<()> {
        let mut buf = Vec::new();

        buf.clear();
        write_position(&mut buf, None)?;
        assert_eq!(buf, b".");

        buf.clear();
        write_position(&mut buf, Some(0))?;
        assert_eq!(buf, b"0");

        Ok(())
    }
}
