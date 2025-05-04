use std::io;

use noodles_bcf as bcf;
use noodles_core::Position;
use noodles_vcf as vcf;

/// A variant record.
#[derive(Clone)]
pub enum Record {
    /// A VCF record.
    Vcf(vcf::Record),
    /// A BCF record.
    Bcf(bcf::Record),
}

impl vcf::variant::Record for Record {
    fn reference_sequence_name<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> io::Result<&'a str> {
        match self {
            Self::Vcf(record) => vcf::variant::Record::reference_sequence_name(record, header),
            Self::Bcf(record) => vcf::variant::Record::reference_sequence_name(record, header),
        }
    }

    fn variant_start(&self) -> Option<io::Result<Position>> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::variant_start(record),
            Record::Bcf(record) => vcf::variant::Record::variant_start(record),
        }
    }

    fn ids(&self) -> Box<dyn vcf::variant::record::Ids + '_> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::ids(record),
            Record::Bcf(record) => vcf::variant::Record::ids(record),
        }
    }

    fn reference_bases(&self) -> Box<dyn vcf::variant::record::ReferenceBases + '_> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::reference_bases(record),
            Record::Bcf(record) => vcf::variant::Record::reference_bases(record),
        }
    }

    fn alternate_bases(&self) -> Box<dyn vcf::variant::record::AlternateBases + '_> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::alternate_bases(record),
            Record::Bcf(record) => vcf::variant::Record::alternate_bases(record),
        }
    }

    fn quality_score(&self) -> Option<io::Result<f32>> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::quality_score(record),
            Record::Bcf(record) => vcf::variant::Record::quality_score(record),
        }
    }

    fn filters(&self) -> Box<dyn vcf::variant::record::Filters + '_> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::filters(record),
            Record::Bcf(record) => vcf::variant::Record::filters(record),
        }
    }

    fn info(&self) -> Box<dyn vcf::variant::record::Info + '_> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::info(record),
            Record::Bcf(record) => vcf::variant::Record::info(record),
        }
    }

    fn samples(&self) -> io::Result<Box<dyn vcf::variant::record::Samples + '_>> {
        match self {
            Record::Vcf(record) => vcf::variant::Record::samples(record),
            Record::Bcf(record) => vcf::variant::Record::samples(record),
        }
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::Vcf(vcf::Record::default())
    }
}
