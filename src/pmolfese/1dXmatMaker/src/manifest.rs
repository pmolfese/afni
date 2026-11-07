// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Stable event-to-column relations, independent of censoring and fitting.
//!
//! One row relates one event to one basis component/column. Ordinary events
//! can share a column; IM events cannot. IDs follow stimulus, run, source-token
//! order, while global IM columns follow AFNI's original token order. These
//! distinct orders are represented explicitly rather than silently sorted.
//!
//! `&mut dyn Write` is a borrowed output interface: a file, buffered writer,
//! or in-memory `Vec<u8>` can receive the same bytes. It resembles a C callback
//! table, but Rust checks that the object stays alive throughout the call.
//! This module neither closes the caller's stream nor owns the Design.
//! The immediately called closure groups std::io errors so we attach one
//! contextual Maker error at the boundary rather than at every write! call.

use crate::{
    design::Design,
    error::{Error, Result},
    model::{ColumnIndex, RunIndex},
};
use std::io::Write;

/// Zero-based stable event ID within this design, not a matrix column number.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EventIndex(pub usize);
/// Zero-based original indexed stimulus, including intervening nuisance inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StimulusIndex(pub usize);
/// Zero-based expanded response index: modulation_set*components+component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BasisIndex(pub usize);

/// One immutable event/basis-to-column mapping, retained even with zero support.
#[derive(Debug)]
pub struct EventMapping {
    /// Stable event identity in stimulus/run/source-token order.
    pub event: EventIndex,
    /// Original CLI stimulus number minus one.
    pub stimulus: StimulusIndex,
    /// Condition identity; initially the same original stimulus identity.
    pub condition: StimulusIndex,
    /// Original run, unaffected by censoring.
    pub run: RunIndex,
    /// Local onset in seconds, before AFNI's float sampling conversion.
    pub onset: f64,
    /// Actual dm duration or fixed stimulus duration, zero for an impulse.
    pub duration: f64,
    /// Expanded response index; AM2 mean components precede modulation sets.
    pub basis: BasisIndex,
    /// Full design column; not a beta-volume index.
    pub column: ColumnIndex,
}

/// Write the version-1 TSV. No item labels are supplied by the M4 CLI, so '-'
/// denotes unspecified item. This table is not Picard's separate trial schema.
pub fn write_manifest(design: &Design, output: &mut dyn Write) -> Result<()> {
    let result = (|| -> std::io::Result<()> {
        writeln!(output, "# 1dXmatMaker events 1")?;
        writeln!(
            output,
            "event\tstimulus\tcondition\trun\tonset\tduration\tbasis\tcolumn\titem"
        )?;
        for m in design.events() {
            writeln!(
                output,
                "{}\t{}\t{}\t{}\t{:.16e}\t{:.16e}\t{}\t{}\t-",
                m.event.0,
                m.stimulus.0,
                m.condition.0,
                m.run.0,
                m.onset,
                m.duration,
                m.basis.0,
                m.column.0
            )?;
        }
        Ok(())
    })();
    result.map_err(|e| Error::resource(format!("cannot write event manifest: {e}")))
}
