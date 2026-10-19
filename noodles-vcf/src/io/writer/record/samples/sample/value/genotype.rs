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

    const VCF_4_0: FileFormat = FileFormat::new(4, 0);

    #[test]
    fn test_write_genotype() -> io::Result<()> {
        let mut buf = Vec::new();

        let genotype: GenotypeBuf = [Allele::new(Some(0), Phasing::Phased)]
            .into_iter()
            .collect();

        buf.clear();
        let header = Header::builder().set_file_format(VCF_4_0).build();
        write_genotype(&mut buf, &header, &&genotype)?;
        assert_eq!(buf, b"0");

        buf.clear();
        let header = Header::builder().set_file_format(VCF_4_4).build();
        write_genotype(&mut buf, &header, &&genotype)?;
        assert_eq!(buf, b"0");

        Ok(())
    }

    #[test]
    fn test_write_alleles_for_vcf_4_0() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            alleles: &[(Option<usize>, Phasing)],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_alleles(buf, VCF_4_0, alleles)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &[(Some(0), Phasing::Phased)], b"0")?;

        buf.clear();
        let alleles = [(Some(0), Phasing::Unphased)];
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_0, &alleles),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        t(
            &mut buf,
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Unphased)],
            b"0/1",
        )?;

        t(
            &mut buf,
            &[(Some(0), Phasing::Phased), (Some(1), Phasing::Phased)],
            b"0|1",
        )?;

        buf.clear();
        let alleles = [(Some(0), Phasing::Unphased), (Some(1), Phasing::Phased)];
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_0, &alleles),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        buf.clear();
        let alleles = [(Some(0), Phasing::Phased), (Some(1), Phasing::Unphased)];
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_0, &alleles),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        t(
            &mut buf,
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Phased),
            ],
            b"0/1|2",
        )?;

        buf.clear();
        let alleles = [
            (Some(0), Phasing::Phased),
            (Some(1), Phasing::Unphased),
            (Some(2), Phasing::Phased),
        ];
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_0, &alleles),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        t(
            &mut buf,
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
            b"./.",
        )?;

        buf.clear();
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_0, &[]),
            Err(e) if e.kind() == io::ErrorKind::InvalidInput
        ));

        Ok(())
    }

    #[test]
    fn test_write_alleles_for_vcf_4_4() -> io::Result<()> {
        fn t(
            buf: &mut Vec<u8>,
            alleles: &[(Option<usize>, Phasing)],
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_alleles(buf, VCF_4_4, alleles)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let mut buf = Vec::new();

        t(&mut buf, &[(Some(0), Phasing::Phased)], b"0")?;
        t(&mut buf, &[(Some(0), Phasing::Unphased)], b"/0")?;
        t(
            &mut buf,
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Unphased)],
            b"0/1",
        )?;
        t(
            &mut buf,
            &[(Some(0), Phasing::Phased), (Some(1), Phasing::Phased)],
            b"0|1",
        )?;
        t(
            &mut buf,
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Phased)],
            b"/0|1",
        )?;
        t(
            &mut buf,
            &[(Some(0), Phasing::Phased), (Some(1), Phasing::Unphased)],
            b"|0/1",
        )?;
        t(
            &mut buf,
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Phased),
            ],
            b"0/1|2",
        )?;
        t(
            &mut buf,
            &[
                (Some(0), Phasing::Phased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Phased),
            ],
            b"|0/1|2",
        )?;
        t(
            &mut buf,
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
            b"./.",
        )?;

        buf.clear();
        assert!(matches!(
            write_alleles(&mut buf, VCF_4_4, &[]),
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
