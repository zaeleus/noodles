//! noodles errors.

use std::fmt;

/// An error kind.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Kind {
    /// An I/O error.
    Io,
    /// An error that does not fall into any other error kind.
    Other,
    /// A parse error.
    Parse,
    /// A conversion error.
    TryFrom,
}

#[derive(Debug)]
enum Context {
    Simple(Kind),
    Custom(Kind, Box<dyn std::error::Error + Send + Sync>),
}

/// A noodles error.
#[derive(Debug)]
pub struct Error {
    context: Context,
}

impl Error {
    /// Creates an error.
    pub fn new<E>(kind: Kind, error: E) -> Self
    where
        E: Into<Box<dyn std::error::Error + Send + Sync>>,
    {
        Self {
            context: Context::Custom(kind, error.into()),
        }
    }

    /// Returns the error kind.
    pub fn kind(&self) -> Kind {
        match &self.context {
            Context::Simple(kind) | Context::Custom(kind, _) => *kind,
        }
    }
}

impl From<Kind> for Error {
    fn from(kind: Kind) -> Self {
        Self {
            context: Context::Simple(kind),
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self::new(Kind::Io, e)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.context {
            Context::Simple(kind) => write!(f, "{:?}", kind),
            Context::Custom(_, e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.context {
            Context::Simple(_) => None,
            Context::Custom(_, e) => Some(&**e),
        }
    }
}
