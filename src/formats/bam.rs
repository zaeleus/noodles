pub use self::cigar::Cigar;
pub use self::data::Data;
pub use self::flag::Flag;
pub use self::quality::Quality;
pub use self::sequence::Sequence;
pub use self::reader::{Reader, Records, References};
pub use self::record::Record;
pub use self::reference::Reference;
pub use self::writer::Writer;

pub mod cigar;
pub mod data;
pub mod flag;
pub mod quality;
pub mod reader;
pub mod record;
pub mod reference;
pub mod sequence;
pub mod writer;

pub static MAGIC_NUMBER: &[u8] = b"BAM\x01";
