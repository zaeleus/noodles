pub use self::cigar::Cigar;
pub use self::data::Data;
pub use self::quality::Quality;
pub use self::sequence::Sequence;
pub use self::reader::{Reader, Records, References};
pub use self::record::{Flag, Record};
pub use self::reference::Reference;

pub mod cigar;
pub mod data;
pub mod quality;
pub mod reader;
pub mod record;
pub mod reference;
pub mod sequence;
