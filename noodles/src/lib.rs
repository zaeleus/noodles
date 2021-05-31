#![deny(missing_docs)]

//! **noodles** is a library for handling various bioinformatics file formats. It currently
//! includes readers and writers for BAM 1.6, BCF 2.2, BGZF, CRAM 3.0, FASTA, FASTQ, GFF3, SAM 1.6,
//! tabix, and VCF 4.3.

#[cfg(feature = "bam")]
pub use noodles_bam as bam;

#[cfg(feature = "bcf")]
pub use noodles_bcf as bcf;

#[cfg(feature = "bgzf")]
pub use noodles_bgzf as bgzf;

#[cfg(feature = "core")]
pub use noodles_core as core;

#[cfg(feature = "cram")]
pub use noodles_cram as cram;

#[cfg(feature = "csi")]
pub use noodles_csi as csi;

#[cfg(feature = "fasta")]
pub use noodles_fasta as fasta;

#[cfg(feature = "fastq")]
pub use noodles_fastq as fastq;

#[cfg(feature = "gff")]
pub use noodles_gff as gff;

#[cfg(feature = "sam")]
pub use noodles_sam as sam;

#[cfg(feature = "tabix")]
pub use noodles_tabix as tabix;

#[cfg(feature = "vcf")]
pub use noodles_vcf as vcf;
