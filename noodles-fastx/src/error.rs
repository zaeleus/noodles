//!

/// Error could be generate by FASTX reading and writing
pub enum Error {
    /// Not the good sequence prefix
    MissingPrefix,
    /// No id
    MissingName,
    /// No sequence
    MissingSequence,
    /// No second description (fastq only)
    MissingSecondDescription,
    /// No or empty quality string (fastq only)
    MissingQuality,
    /// Partial Record
    PartialRecord,
}

// Convert an std::io::Error in Error is impossible this clippy warning isn't intrusting here
#[allow(clippy::from_over_into)]
impl Into<std::io::Error> for Error {
    fn into(self) -> std::io::Error {
        match self {
            Error::MissingPrefix => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The prefix ('>' or '@') is missing",
            ),
            Error::MissingSequence => {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "The sequence is missing")
            }
            Error::MissingName => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The sequence name is missing",
            ),
            Error::MissingSecondDescription => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The second description is missing",
            ),
            Error::MissingQuality => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The quality string is missing",
            ),
            Error::PartialRecord => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The last fastq record isn't complete",
            ),
        }
    }
}
