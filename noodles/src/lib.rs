#![warn(missing_docs)]

//! **noodles** is a library for handling various bioinformatics file formats. It currently
//! includes readers and writers for BAM 1.6, BCF 2.2, BED, BGZF, CRAM 3.0, CSI, FASTA, FASTQ,
//! GFF3, GTF 2.2, SAM 1.6, tabix, and VCF 4.3.

#[cfg(feature = "bam")]
#[doc(inline)]
pub use noodles_bam as bam;

#[cfg(feature = "bcf")]
#[doc(inline)]
pub use noodles_bcf as bcf;

#[cfg(feature = "bed")]
#[doc(inline)]
pub use noodles_bed as bed;

#[cfg(feature = "bgzf")]
#[doc(inline)]
pub use noodles_bgzf as bgzf;

#[cfg(feature = "core")]
#[doc(inline)]
pub use noodles_core as core;

#[cfg(feature = "cram")]
#[doc(inline)]
pub use noodles_cram as cram;

#[cfg(feature = "csi")]
#[doc(inline)]
pub use noodles_csi as csi;

#[cfg(feature = "fasta")]
#[doc(inline)]
pub use noodles_fasta as fasta;

#[cfg(feature = "fastq")]
#[doc(inline)]
pub use noodles_fastq as fastq;

#[cfg(feature = "gff")]
#[doc(inline)]
pub use noodles_gff as gff;

#[cfg(feature = "gtf")]
#[doc(inline)]
pub use noodles_gtf as gtf;

#[cfg(feature = "sam")]
#[doc(inline)]
pub use noodles_sam as sam;

#[cfg(feature = "tabix")]
#[doc(inline)]
pub use noodles_tabix as tabix;

#[cfg(feature = "vcf")]
#[doc(inline)]
pub use noodles_vcf as vcf;
