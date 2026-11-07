// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Stage complete outputs, then publish with atomic no-replace hard links.
//!
//! A plain rename can overwrite an existing destination on Unix. A hard link
//! on the same filesystem instead fails when the destination exists. Each
//! Stage owns its temporary file and remembers whether it published a link;
//! Drop cleans up on any early return. Two files cannot be atomically committed
//! as one filesystem transaction, but ordinary failures trigger rollback.
//!
//! Supported release platforms are Unix (Linux/macOS). No unsafe code or FFI
//! is used. Rollback checks inode identity so it will not remove a replacement
//! created by another process. As with ordinary filesystem programs, hostile
//! concurrent directory replacement is not a supported transaction model.
//!
//! Newcomer reading order: destinations validates all requested names;
//! create_stage obtains exclusive temporary files; publish writes every stream
//! and then publishes links. Stage is the cleanup guard that connects them.
//! Its Drop method runs on ordinary scope exit AND early returns using `?`.
//! The atomic counter only chooses temporary names; it does not make a
//! multi-file transaction atomic or synchronize scientific computations.

use crate::error::{Error, Result};
use std::{
    collections::HashSet,
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    os::unix::fs::MetadataExt,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

static SERIAL: AtomicU64 = AtomicU64::new(0);

struct Stage {
    temporary: PathBuf,
    destination: PathBuf,
    identity: (u64, u64),
    published: bool,
    committed: bool,
}

impl Drop for Stage {
    fn drop(&mut self) {
        if self.published && !self.committed {
            // Never follow a symlink during rollback, and never remove a new
            // inode that has replaced our link since publication.
            if let Ok(metadata) = fs::symlink_metadata(&self.destination)
                && metadata.is_file()
                && (metadata.dev(), metadata.ino()) == self.identity
            {
                let _ = fs::remove_file(&self.destination);
            }
        }
        let _ = fs::remove_file(&self.temporary);
    }
}

/// Resolve parent aliases and reject every existing output, including dangling
/// symlinks. Canonicalizing only the parent works when the output is new.
fn destinations(paths: &[&Path]) -> Result<Vec<PathBuf>> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    for path in paths {
        let name = path
            .file_name()
            .filter(|n| *n != "-")
            .ok_or_else(|| Error::input(format!("invalid output filename '{}'", path.display())))?;
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let canonical = fs::canonicalize(parent)
            .map_err(|e| Error::io("cannot resolve output directory", parent, e))?;
        let destination = canonical.join(name);
        if !seen.insert(destination.clone()) {
            return Err(Error::input("output paths resolve to the same destination"));
        }
        match fs::symlink_metadata(&destination) {
            Ok(_) => {
                return Err(Error::input(format!(
                    "output already exists: {}",
                    destination.display()
                )));
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(Error::io("cannot inspect output", &destination, e)),
        }
        result.push(destination);
    }
    Ok(result)
}

/// Create a private staging file beside the output using exclusive creation.
/// Predictable names cannot cause overwrites: create_new rejects even symlinks.
fn create_stage(destination: PathBuf) -> Result<(Stage, File)> {
    let parent = destination
        .parent()
        .ok_or_else(|| Error::input("output has no parent"))?;
    for _ in 0..100 {
        let serial = SERIAL.fetch_add(1, Ordering::Relaxed);
        let temporary = parent.join(format!(".1dXmatMaker-{}-{serial}.tmp", std::process::id()));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => {
                let metadata = match file.metadata() {
                    Ok(metadata) => metadata,
                    Err(e) => {
                        let _ = fs::remove_file(&temporary);
                        return Err(Error::io("cannot inspect staging file", &temporary, e));
                    }
                };
                return Ok((
                    Stage {
                        temporary,
                        destination,
                        identity: (metadata.dev(), metadata.ino()),
                        published: false,
                        committed: false,
                    },
                    file,
                ));
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(Error::io("cannot stage output", &temporary, e)),
        }
    }
    Err(Error::resource(
        "could not obtain an unused staging filename",
    ))
}

/// Write all files completely before linking any destination into place.
/// `write_file` receives an output index and a borrowed writer, not ownership
/// of the file. A writer error returns immediately and drops every stage.
pub(crate) fn publish(
    paths: &[&Path],
    mut write_file: impl FnMut(usize, &mut dyn Write) -> Result<()>,
) -> Result<()> {
    let destinations = destinations(paths)?;
    let mut stages = Vec::new();
    for (index, destination) in destinations.into_iter().enumerate() {
        let (stage, file) = create_stage(destination)?;
        let mut writer = BufWriter::new(file);
        write_file(index, &mut writer)?;
        writer
            .flush()
            .map_err(|e| Error::io("cannot flush staging file", &stage.temporary, e))?;
        writer
            .get_ref()
            .sync_all()
            .map_err(|e| Error::io("cannot sync staging file", &stage.temporary, e))?;
        stages.push(stage);
    }
    for stage in &mut stages {
        fs::hard_link(&stage.temporary, &stage.destination).map_err(|e| {
            if e.kind() == std::io::ErrorKind::AlreadyExists {
                Error::input(format!(
                    "output appeared during publication: {}",
                    stage.destination.display()
                ))
            } else {
                Error::io("cannot publish output", &stage.destination, e)
            }
        })?;
        stage.published = true;
    }
    for stage in &mut stages {
        stage.committed = true;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    //! Fault injection tests the publication transaction, not numeric formatting.
    use super::*;

    struct Scratch(PathBuf);
    impl Scratch {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "xmat-publication-{}-{}",
                std::process::id(),
                SERIAL.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).expect("test scratch directory");
            Self(path)
        }
    }
    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn competing_second_output_rolls_back_only_our_first_file() {
        let scratch = Scratch::new();
        let first = scratch.0.join("first");
        let second = scratch.0.join("second");
        let result = publish(&[&first, &second], |index, writer| {
            writer.write_all(b"our matrix").unwrap();
            // This appears AFTER preflight but BEFORE publication. The first
            // link succeeds; the second must fail without replacing this file.
            if index == 1 {
                fs::write(&second, b"another writer").unwrap();
            }
            Ok(())
        });
        assert!(result.is_err());
        assert!(!first.exists());
        assert_eq!(fs::read(&second).unwrap(), b"another writer");
        assert_eq!(fs::read_dir(&scratch.0).unwrap().count(), 1);
    }

    #[test]
    fn manifest_collision_rolls_back_both_matrices() {
        // M2 adds a third output. A competitor can create its destination
        // after preflight; rollback must then remove BOTH matrix links while
        // leaving that competitor's manifest and every preexisting file alone.
        let scratch = Scratch::new();
        let retained = scratch.0.join("retained");
        let full = scratch.0.join("full");
        let manifest = scratch.0.join("events");
        let result = publish(&[&retained, &full, &manifest], |index, writer| {
            writer.write_all(b"complete staged content").unwrap();
            if index == 2 {
                fs::write(&manifest, b"competing manifest").unwrap();
            }
            Ok(())
        });
        assert!(result.is_err());
        assert!(!retained.exists() && !full.exists());
        assert_eq!(fs::read(&manifest).unwrap(), b"competing manifest");
        assert_eq!(fs::read_dir(&scratch.0).unwrap().count(), 1);
    }

    #[test]
    fn write_failure_cleans_all_stages_before_publication() {
        let scratch = Scratch::new();
        let first = scratch.0.join("first");
        let second = scratch.0.join("second");
        let result = publish(&[&first, &second], |index, writer| {
            writer.write_all(b"partial data").unwrap();
            if index == 1 {
                return Err(Error::resource("injected disk failure"));
            }
            Ok(())
        });
        assert!(result.is_err());
        assert_eq!(fs::read_dir(&scratch.0).unwrap().count(), 0);
    }
}
