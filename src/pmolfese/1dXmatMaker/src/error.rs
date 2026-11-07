// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Errors cross library boundaries as values, not process exits or panics.
//!
//! Rust's `Result<T, Error>` contains either a successful `T` or an `Error`.
//! Callers use `?` to return an error early while owned resources are cleaned
//! up automatically. Only main.rs decides the process exit code.
//!
//! The `Into<String>` parameter accepts either borrowed literal text or an
//! owned String and gives Error its own message. `map_err` translates one
//! error type to another; it does not suppress failure. The Display trait is
//! the formatting contract used by main's eprintln!, while std::error::Error
//! lets other Rust applications handle these errors using standard interfaces.

use std::{fmt, io, path::Path};

/// Classification used by the CLI; library users can inspect it directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    /// Invalid or unsupported input, including an output collision (exit 2).
    Input,
    /// File-system, allocation, or internal resource failure (exit 1).
    Resource,
}

/// An error with human-readable file/option/dimension context.
#[derive(Debug)]
pub struct Error {
    kind: ErrorKind,
    message: String,
}

impl Error {
    /// Describe invalid data without terminating the calling process.
    pub fn input(message: impl Into<String>) -> Self {
        Self {
            kind: ErrorKind::Input,
            message: message.into(),
        }
    }

    /// Describe a failure to obtain storage or complete an operation.
    pub fn resource(message: impl Into<String>) -> Self {
        Self {
            kind: ErrorKind::Resource,
            message: message.into(),
        }
    }

    /// Add the operation and path to an operating-system error.
    pub fn io(action: &str, path: &Path, cause: io::Error) -> Self {
        Self::resource(format!("{action} '{}': {cause}", path.display()))
    }

    /// Return the error category without borrowing its message.
    pub fn kind(&self) -> ErrorKind {
        self.kind
    }

    /// Return the documented CLI failure code, 2 for input and 1 otherwise.
    pub fn exit_code(&self) -> u8 {
        match self.kind {
            ErrorKind::Input => 2,
            ErrorKind::Resource => 1,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.message)
    }
}
impl std::error::Error for Error {}

/// Library result type; `Ok(value)` is success and `Err(error)` is failure.
pub type Result<T> = std::result::Result<T, Error>;

/// Reserve space before extending a scientific buffer; report capacity failure.
///
/// This helper avoids infallible large `vec![...; user_size]` allocations.
/// The OS can still kill a process under system-wide memory pressure.
pub(crate) fn reserve<T>(values: &mut Vec<T>, additional: usize) -> Result<()> {
    values
        .try_reserve_exact(additional)
        .map_err(|e| Error::resource(format!("cannot allocate {additional} elements: {e}")))
}
