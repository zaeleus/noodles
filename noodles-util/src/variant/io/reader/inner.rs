use std::io::{self, BufRead};

use noodles_bcf as bcf;
use noodles_vcf::{self as vcf, variant::io::Read};

use crate::variant::Record;

pub(super) enum Inner<R> {
    Vcf(vcf::io::Reader<R>),
    Bcf(bcf::io::Reader<R>),
}

impl<R> Inner<R>
where
    R: BufRead,
{
    pub(super) fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Vcf(reader) => reader.read_header(),
            Self::Bcf(reader) => reader.read_header(),
        }
    }

    pub(super) fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        match self {
            Inner::Vcf(reader) => {
                if !matches!(record, Record::Vcf(_)) {
                    *record = Record::Vcf(vcf::Record::default());
                }

                if let Record::Vcf(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
            Inner::Bcf(reader) => {
                if !matches!(record, Record::Bcf(_)) {
                    *record = Record::Bcf(bcf::Record::default());
                }

                if let Record::Bcf(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
        }
    }

    pub(super) fn records<'a>(
        &'a mut self,
        header: &'a vcf::Header,
    ) -> impl Iterator<Item = io::Result<Box<dyn vcf::variant::Record>>> + 'a {
        match self {
            Inner::Vcf(reader) => reader.variant_records(header),
            Inner::Bcf(reader) => reader.variant_records(header),
        }
    }
}
