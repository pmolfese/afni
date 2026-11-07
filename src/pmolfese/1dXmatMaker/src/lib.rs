// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Independent AFNI matrix generation: acquisition -> full design -> row map.
//!
//! Start with [`model::Acquisition`], then [`design::build_design`], then
//! [`xmat::write_outputs`]. The CLI is a separate adapter; these modules neither
//! invoke AFNI nor inspect voxel data. The response catalog, EXPR, timing,
//! modulation, sampled lags and GLTs form the ordinary path. For versioned
//! sampled-HRF collections, start with [`candidates`] and [`sampled`].
//!
//! A `Design` owns one column-major full-timeline matrix. Censoring is an index
//! map into it, so writing retained and full outputs cannot accidentally build
//! two different versions of the same regressor. All values are finite f64.
//!
//! Build a two-run baseline without files or the CLI:
//! ```
//! use afni_xmat_maker::{design::build_design, model::Acquisition};
//! # fn main() -> afni_xmat_maker::error::Result<()> {
//! let acquisition = Acquisition::new(vec![20, 28], 1.0, 1, None)?;
//! let design = build_design(acquisition, Vec::new(), None)?;
//! assert_eq!(design.columns().len(), 4); // constant + slope for each run
//! assert_eq!(design.retained_rows().len(), 48);
//! # Ok(())
//! # }
//! ```
//!
//! This package has no afni-io or afni-core dependency. Similar names in
//! source comments describe scientific/format provenance, not external calls.
//! If either crate is adopted, call-site documentation must state why it is
//! needed, the input/output shape and units, who owns the buffers, how errors
//! are translated, and which compatibility assumptions the call introduces.

//! Build two individual events and inspect their stable mapping:
//! ```
//! use afni_xmat_maker::{basis::Basis, design::{build_with_events, Source,
//!     EventStimulus, EventMode}, model::Acquisition, timing::{Timing, TimingMode}};
//! # fn main() -> afni_xmat_maker::error::Result<()> {
//! let acquisition = Acquisition::new(vec![40], 1.0, 0, None)?;
//! let times = Timing::parse("2.5 18", "example", TimingMode::Local, &acquisition)?;
//! let trials = EventStimulus::new("trial", Basis::parse("GAM")?, times,
//!     EventMode::Individual)?;
//! let design = build_with_events(acquisition, vec![Source::Events(trials)], None)?;
//! assert_eq!(design.events()[0].column.0, 1); // column 0 is the intercept
//! assert_eq!(design.events()[1].event.0, 1); // each IM event keeps its own ID
//! # Ok(())
//! # }
//! ```

pub mod basis;
mod basis_extended;
pub mod candidates;
pub mod censor;
pub mod contrast;
pub mod design;
pub mod error;
mod expression;
mod expression_math;
pub mod input;
pub mod manifest;
pub mod model;
mod publication;
pub mod sampled;
pub mod table;
pub mod timing;
pub mod xmat;

/// Version of the independent generator, distinct from AFNI's version.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Pinned scientific behavior reference; provenance is in
/// afni/compatibility/reference/manifest.json inside fixtures.tgz.
pub const REFERENCE_REVISION: &str = "2664fdd5de7425f658c6d901765fbb4c291693fa";
