bitflags::bitflags! {
    /// SAM record fields.
    #[derive(Default)]
    pub struct Fields: u16 {
        /// Read name (`QNAME`).
        const READ_NAME = 0x01;
        /// Flags (`FLAG`).
        const FLAGS = 0x02;
        /// Reference sequecne ID (`RNAME` equiv.).
        const REFERENCE_SEQUENCE_ID = 0x04;
        /// Position (`POS`).
        const ALIGNMENT_START = 0x08;
        /// Mapping quality (`MAPQ`).
        const MAPPING_QUALITY = 0x10;
        /// CIGAR operations (`CIGAR`).
        const CIGAR = 0x20;
        /// Mate reference sequence ID (`RNEXT` equiv.).
        const MATE_REFERENCE_SEQUENCE_ID = 0x40;
        /// Mate position (`PNEXT`).
        const MATE_ALIGNMENT_START = 0x80;
        /// Template length (`TLEN`).
        const TEMPLATE_LENGTH = 0x0100;
        /// Sequence (`SEQ`).
        const SEQUENCE = 0x0200;
        /// Quality scores (`QUAL`).
        const QUALITY_SCORES = 0x0400;
        /// Data.
        const DATA = 0x0800;
    }
}
