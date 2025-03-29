//! **noodles-util** are utilities for working with noodles. Currently, this consists of a unified
//! interface for reading and writing alignment (BAM/CRAM/SAM) and variant (VCF/BCF) data.

#[cfg(feature = "alignment")]
pub mod alignment;

#[cfg(feature = "variant")]
pub mod variant;
