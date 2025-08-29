use std::io::{self, BufReader, Read};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, variant::io::Read as _};

use crate::variant::Record;

pub(super) enum Inner<R> {
    Bcf(bcf::io::Reader<bgzf::io::Reader<BufReader<R>>>),
    BcfRaw(bcf::io::Reader<BufReader<R>>),
    Vcf(vcf::io::Reader<BufReader<R>>),
    VcfGz(vcf::io::Reader<bgzf::io::Reader<BufReader<R>>>),
}

impl<R> Inner<R>
where
    R: Read,
{
    pub(super) fn read_header(&mut self) -> io::Result<vcf::Header> {
        match self {
            Self::Bcf(reader) => reader.read_header(),
            Self::BcfRaw(reader) => reader.read_header(),
            Self::Vcf(reader) => reader.read_header(),
            Self::VcfGz(reader) => reader.read_header(),
        }
    }

    pub(super) fn read_record(&mut self, record: &mut Record) -> io::Result<usize> {
        match self {
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
            Inner::BcfRaw(reader) => {
                if !matches!(record, Record::Bcf(_)) {
                    *record = Record::Bcf(bcf::Record::default());
                }

                if let Record::Bcf(r) = record {
                    reader.read_record(r)
                } else {
                    unreachable!();
                }
            }
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
            Inner::VcfGz(reader) => {
                if !matches!(record, Record::Vcf(_)) {
                    *record = Record::Vcf(vcf::Record::default());
                }

                if let Record::Vcf(r) = record {
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
            Inner::Bcf(reader) => reader.variant_records(header),
            Inner::BcfRaw(reader) => reader.variant_records(header),
            Inner::Vcf(reader) => reader.variant_records(header),
            Inner::VcfGz(reader) => reader.variant_records(header),
        }
    }
}
