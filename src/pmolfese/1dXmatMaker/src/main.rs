// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Thin process adapter: parse arguments, call the library, report exit status.
//!
//! `ExitCode` keeps normal errors out of the panic mechanism. In particular,
//! returning from main lets Rust drop owned buffers and staging guards; calling
//! process::exit deep inside a numerical function would bypass that cleanup.
//!
//! The binary and library are separate Cargo targets. `mod cli` compiles this
//! binary's argument adapter; the adapter imports the reusable library by its
//! package name. A library caller can therefore build matrices without parsing
//! argv, printing messages, or launching 1dXmatMaker as a subprocess.

mod candidate_cli;
mod cli;
use std::process::ExitCode;

fn main() -> ExitCode {
    match cli::run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("** 1dXmatMaker: {error}");
            ExitCode::from(error.exit_code())
        }
    }
}
