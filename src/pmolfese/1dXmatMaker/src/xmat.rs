// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Serialize a validated design as an ordinary AFNI NIML-comment matrix.
//!
//! The writer borrows the design; it cannot change values, labels, or row maps.
//! Numeric rows are streamed from column-major storage without allocating a
//! second censored matrix. See afni/compatibility/reference inside fixtures.tgz
//! for the independent AFNI writer examples and tests for actual AFNI reader
//! interoperability.
//!
//! Format compatibility is not a library call: no afni-io NIML writer or
//! afni-core object is used here. `attribute` escapes text into the established
//! AFNI header grammar, and the shared C reader is exercised by tests only.
//! `write_matrix` accepts any Write stream, including a test memory buffer;
//! `write_outputs` adds filesystem staging and owns publication responsibilities.
//! This split keeps a formatting routine from deciding output filenames or
//! replacing existing files as a side effect of serializing a Design.

use crate::{
    REFERENCE_REVISION, VERSION,
    design::Design,
    error::{Error, Result},
    manifest::write_manifest,
    model::{ColumnIndex, FullRow},
    publication,
};
use std::{io::Write, path::Path};

/// Which rows to serialize from the same immutable full design.
#[derive(Debug, Clone, Copy)]
pub enum RowSelection {
    /// GoodList selects the retained original row indices.
    Retained,
    /// Identity GoodList writes all original rows, including censored rows.
    Full,
    /// Full rows with one additional unit spike per excluded original row.
    Regcensored,
}

/// Write one xmat stream; partial writes on an arbitrary writer are possible.
///
/// Use write_outputs for no-overwrite staged file publication. command_line
/// is provenance text only: it is escaped, never executed. Signal-only stimulus
/// ranges and original-index basis metadata are emitted when events exist.
pub fn write_matrix(
    design: &Design,
    selection: RowSelection,
    command_line: &str,
    output: &mut dyn Write,
) -> Result<()> {
    let acquisition = design.acquisition();
    // Spikes use original row IDs. Retained rows stay untouched; each excluded
    // row gets its own nuisance coefficient in the full-row companion.
    let excluded: Vec<_> = if matches!(selection, RowSelection::Regcensored) {
        (0..acquisition.nrows())
            .filter(|r| {
                design
                    .retained_rows()
                    .binary_search_by_key(r, |row| row.0)
                    .is_err()
            })
            .collect()
    } else {
        Vec::new()
    };
    let original_columns = design.columns().len();
    let ncols = original_columns + excluded.len();
    let count = match selection {
        RowSelection::Full | RowSelection::Regcensored => acquisition.nrows(),
        RowSelection::Retained => design.retained_rows().len(),
    };
    // Keep serialization errors contextual without exposing I/O internals in
    // numerical APIs. File setup, flush, and publication failures separately
    // include their path in the publication layer.
    let io_result = (|| -> std::io::Result<()> {
        writeln!(output, "# <matrix")?;
        attribute(output, "ni_type", &format!("{ncols}*double"))?;
        attribute(output, "ni_dimen", &count.to_string())?;
        let mut labels = design
            .columns()
            .iter()
            .map(|c| c.label())
            .collect::<Vec<_>>()
            .join(" ; ");
        for row in &excluded {
            labels.push_str(&format!(" ; cens#{row}"));
        }
        attribute(output, "ColumnLabels", &labels)?;
        let mut groups = design
            .columns()
            .iter()
            .map(|c| c.group().to_string())
            .collect::<Vec<_>>()
            .join(",");
        for _ in &excluded {
            groups.push_str(",0");
        }
        attribute(output, "ColumnGroups", &groups)?;
        attribute(output, "RowTR", &acquisition.tr().to_string())?;
        // Lists need not use AFNI's compressed syntax. Streaming avoids a
        // potentially large temporary string for the full identity GoodList.
        write!(output, "#  GoodList = \"")?;
        for position in 0..count {
            if position > 0 {
                write!(output, ",")?;
            }
            write!(output, "{}", selected_row(design, selection, position).0)?;
        }
        writeln!(output, "\"")?;
        attribute(output, "NRowFull", &acquisition.nrows().to_string())?;
        let starts = acquisition
            .run_starts()
            .iter()
            .map(|r| r.0.to_string())
            .collect::<Vec<_>>()
            .join(",");
        attribute(output, "RunStart", &starts)?;
        if !design.stimuli().is_empty() {
            let stimuli = design.stimuli();
            attribute(output, "Nstim", &stimuli.len().to_string())?;
            attribute(
                output,
                "StimBots",
                &stimuli
                    .iter()
                    .map(|s| s.first.0.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            )?;
            attribute(
                output,
                "StimTops",
                &stimuli
                    .iter()
                    .map(|s| s.last.0.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
            )?;
            attribute(
                output,
                "StimLabels",
                &stimuli
                    .iter()
                    .map(|s| s.label.as_str())
                    .collect::<Vec<_>>()
                    .join(" ; "),
            )?;
            // BasisNstim includes nuisance slots and appended ortvec columns;
            // the six-digit suffix is the ORIGINAL 1-based stimulus number.
            if stimuli.iter().any(|s| !s.formula.is_empty()) {
                attribute(output, "BasisNstim", &design.basis_nstim().to_string())?;
            }
            for stimulus in stimuli.iter().filter(|s| !s.formula.is_empty()) {
                let suffix = format!("{:06}", stimulus.index.0 + 1);
                attribute(
                    output,
                    &format!("BasisOption_{suffix}"),
                    stimulus.mode.option(),
                )?;
                attribute(output, &format!("BasisName_{suffix}"), &stimulus.basis_name)?;
                attribute(output, &format!("BasisFormula_{suffix}"), &stimulus.formula)?;
                attribute(
                    output,
                    &format!("BasisColumns_{suffix}"),
                    &format!("{}:{}", stimulus.first.0, stimulus.last.0),
                )?;
            }
        }
        // AFNI associates GLTs only with the primary retained matrix. Full and
        // spike-augmented companions deliberately omit these coefficient maps.
        if matches!(selection, RowSelection::Retained) && !design.contrasts.is_empty() {
            attribute(output, "Nglt", &design.contrasts.len().to_string())?;
            attribute(
                output,
                "GltLabels",
                &design
                    .contrasts
                    .iter()
                    .map(|g| g.label.as_str())
                    .collect::<Vec<_>>()
                    .join(" ; "),
            )?;
            for (index, glt) in design.contrasts.iter().enumerate() {
                let mut values = vec![glt.rows.len().to_string(), ncols.to_string()];
                values.extend(glt.rows.iter().flatten().map(|v| v.to_string()));
                attribute(output, &format!("GltMatrix_{index:06}"), &values.join(","))?;
            }
        }
        attribute(output, "CommandLine", command_line)?;
        attribute(output, "XmatMakerVersion", VERSION)?;
        attribute(output, "XmatMakerReference", REFERENCE_REVISION)?;
        writeln!(output, "# >")?;
        for position in 0..count {
            let row = selected_row(design, selection, position);
            for column in 0..ncols {
                // Row/column bounds were validated by the builder and loops.
                let value = if column < original_columns {
                    design
                        .value(row, ColumnIndex(column))
                        .ok_or_else(|| std::io::Error::other("internal matrix shape mismatch"))?
                } else {
                    f64::from(row.0 == excluded[column - original_columns])
                };
                if column > 0 {
                    write!(output, " ")?;
                }
                // 17 significant digits round-trip binary64. Scientific
                // notation keeps very small/large values compact and finite.
                write!(output, "{value:.16e}")?;
            }
            writeln!(output)?;
        }
        writeln!(output, "# </matrix>")
    })();
    io_result.map_err(|e| Error::resource(format!("cannot write xmat stream: {e}")))
}

/// Map a valid output position to the corresponding original row.
fn selected_row(design: &Design, selection: RowSelection, position: usize) -> FullRow {
    match selection {
        RowSelection::Full | RowSelection::Regcensored => FullRow(position),
        RowSelection::Retained => design.retained_rows()[position],
    }
}

/// Escape every attribute value, including command provenance and labels.
/// Numeric values use this same path so no header can inject a new attribute.
fn attribute(output: &mut dyn Write, name: &str, value: &str) -> std::io::Result<()> {
    write!(output, "#  {name} = \"")?;
    for character in value.chars() {
        match character {
            '&' => write!(output, "&amp;")?,
            '<' => write!(output, "&lt;")?,
            '>' => write!(output, "&gt;")?,
            '"' => write!(output, "&quot;")?,
            '\'' => write!(output, "&apos;")?,
            // AFNI NIML does not promise all XML numeric character entities.
            // Flatten command-line control characters for a single-line header.
            c if c.is_control() => write!(output, " ")?,
            c => write!(output, "{c}")?,
        }
    }
    writeln!(output, "\"")
}

/// Publish retained output and an optional full companion without overwrites.
///
/// Parent directories must exist. All outputs are staged before publication;
/// failure rolls back only links created by this invocation. See publication.rs
/// for filesystem transaction limits. Omitting full writes exactly one file.
pub fn write_outputs(
    design: &Design,
    retained: &Path,
    full: Option<&Path>,
    command_line: &str,
) -> Result<()> {
    write_outputs_with_manifest(design, retained, full, None, command_line)
}

/// Publish the retained matrix, optional full matrix, and optional event TSV
/// through one staging/rollback operation. No default sidecar is created.
/// Requesting a manifest without any event source is an input error.
pub fn write_outputs_with_manifest(
    design: &Design,
    retained: &Path,
    full: Option<&Path>,
    manifest: Option<&Path>,
    command_line: &str,
) -> Result<()> {
    write_all_outputs(design, retained, full, manifest, None, command_line)
}

/// Stage every requested artifact before publishing any path. Keeping the
/// (path,kind) pairs together prevents optional companions from shifting an
/// index onto the wrong serializer. No output is overwritten on failure.
pub fn write_all_outputs(
    design: &Design,
    retained: &Path,
    full: Option<&Path>,
    manifest: Option<&Path>,
    regcen: Option<&Path>,
    command_line: &str,
) -> Result<()> {
    if manifest.is_some() && design.stimuli().iter().all(|s| s.formula.is_empty()) {
        return Err(Error::input("-event_manifest requires an event stimulus"));
    }
    let mut artifacts = vec![(retained, Some(RowSelection::Retained))];
    if let Some(p) = full {
        artifacts.push((p, Some(RowSelection::Full)));
    }
    if let Some(p) = regcen {
        artifacts.push((p, Some(RowSelection::Regcensored)));
    }
    if let Some(p) = manifest {
        artifacts.push((p, None));
    }
    let paths = artifacts.iter().map(|(p, _)| *p).collect::<Vec<_>>();
    publication::publish(&paths, |index, writer| match artifacts[index].1 {
        Some(selection) => write_matrix(design, selection, command_line, writer),
        None => write_manifest(design, writer),
    })
}
