//! **noodles-gtf** handles the reading and writing of the Gene Transfer Format (GTF).

pub mod io;
pub mod line;
mod line_buf;
pub mod record;

pub use self::{line::Line, line_buf::LineBuf, record::Record};

#[deprecated(since = "0.32.0", note = "Use `noodles_gtf::io::Reader` instead.")]
pub use self::io::Reader;

#[deprecated(since = "0.32.0", note = "Use `noodles_gtf::io::Writer` instead.")]
pub use self::io::Writer;
