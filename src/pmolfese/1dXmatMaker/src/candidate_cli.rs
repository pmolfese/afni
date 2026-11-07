// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! CLI adapter for aligned sampled-HRF collections. The ordinary AFNI-style
//! parser remains in cli.rs; this mode has explicit events and a library so
//! it cannot accidentally mix a candidate campaign with unrelated stim_times.
//! Parsing gathers owned settings first. The local library validates waveforms
//! and identities, then publishes all outputs with rollback on ordinary errors.
use afni_xmat_maker::{
    candidates::{self, Request},
    design::{Nuisance, NuisanceKind},
    error::{Error, Result},
    model::Acquisition,
    table::{read_censor, read_table},
};
use std::{collections::HashMap, path::PathBuf};

/// Consume the candidate-mode option list. Duplicate singleton or unknown
/// options fail; ortvec alone may repeat. No AFNI environment is consulted.
pub(crate) fn run(args: &[String], command: &str) -> Result<()> {
    let mut settings = HashMap::new();
    let mut ortvec = Vec::new();
    let mut cursor = 0;
    while cursor < args.len() {
        let option = &args[cursor];
        cursor += 1;
        let value = args
            .get(cursor)
            .ok_or_else(|| Error::input(format!("{option} needs a value")))?;
        cursor += 1;
        if option == "-ortvec" {
            let label = args
                .get(cursor)
                .ok_or_else(|| Error::input("ortvec needs a label"))?;
            cursor += 1;
            ortvec.push((PathBuf::from(value), label.clone()));
            continue;
        }
        if !matches!(
            option.as_str(),
            "-candidate_library"
                | "-events"
                | "-candidate_prefix"
                | "-run_lengths"
                | "-TR"
                | "-polort"
                | "-censor"
                | "-candidate_ids"
                | "-hrf_normalization"
        ) {
            return Err(Error::input(format!(
                "unsupported candidate option '{option}'"
            )));
        }
        if settings.insert(option.as_str(), value.as_str()).is_some() {
            return Err(Error::input(format!("duplicate {option}")));
        }
    }
    let required = |key| {
        settings
            .get(key)
            .copied()
            .ok_or_else(|| Error::input(format!("candidate mode requires {key}")))
    };
    let lengths = required("-run_lengths")?
        .split(',')
        .map(|s| {
            s.parse::<usize>()
                .map_err(|_| Error::input("invalid run length"))
        })
        .collect::<Result<Vec<_>>>()?;
    let tr = required("-TR")?
        .parse::<f64>()
        .map_err(|_| Error::input("invalid TR"))?;
    let polort = required("-polort")?
        .parse::<i32>()
        .map_err(|_| Error::input("candidate polort must be explicit integer"))?;
    let acquisition = Acquisition::new(lengths, tr, polort, None)?;
    let trials = candidates::load_trials(&PathBuf::from(required("-events")?), &acquisition)?;
    let mut library = candidates::load_library(
        required("-candidate_library")?,
        trials[0].duration,
        tr,
        settings
            .get("-hrf_normalization")
            .copied()
            .unwrap_or("peak"),
    )?;
    if let Some(ids) = settings.get("-candidate_ids") {
        let requested = ids
            .split(',')
            .map(|v| {
                v.parse::<usize>()
                    .map_err(|_| Error::input("candidate_ids needs comma-separated stable IDs"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut seen = std::collections::HashSet::new();
        if requested
            .iter()
            .any(|id| !seen.insert(*id) || !library.iter().any(|c| c.id == *id))
        {
            return Err(Error::input("unknown or repeated candidate ID"));
        }
        library.retain(|c| seen.contains(&c.id));
    }
    let nuisance = ortvec
        .into_iter()
        .map(|(path, label)| {
            Nuisance::new(
                &label,
                NuisanceKind::Ortvec,
                read_table(&path, acquisition.nrows())?,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let keep = settings
        .get("-censor")
        .map(|p| read_censor(&PathBuf::from(p), acquisition.nrows()))
        .transpose()?;
    let count = library.len();
    let prefix = PathBuf::from(required("-candidate_prefix")?);
    candidates::generate(
        Request {
            acquisition,
            trials,
            library,
            nuisance,
            keep,
            prefix: prefix.clone(),
        },
        command,
    )?;
    eprintln!(
        "++ wrote {count} aligned trial/ONOFF candidates at {}",
        prefix.display()
    );
    Ok(())
}
