//! Editable presentation plus a replayable, fingerprinted source recipe.
//! Computed samples and render caches are deliberately not serialized.
#![cfg_attr(not(feature = "gui"), allow(dead_code))]
use crate::{
    cli,
    model::{self, Bounds, Figure},
};
use anyhow::{Context, Result, ensure};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    io::{Read, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Panel {
    pub id: usize,
    pub title: String,
    pub caption: String,
    pub xlabel: String,
    pub ylabel: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Presentation {
    pub title: String,
    #[serde(default)]
    pub subtitle: Option<String>,
    pub caption: String,
    pub columns: usize,
    pub panels: Vec<Panel>,
}
impl Presentation {
    pub fn capture(fig: &Figure) -> Self {
        Self {
            title: fig.title.clone(),
            subtitle: Some(fig.subtitle.clone()),
            caption: fig.caption.clone(),
            columns: fig.columns,
            panels: fig
                .plots
                .iter()
                .map(|p| Panel {
                    id: p.id,
                    title: p.title.clone(),
                    caption: p.caption.clone(),
                    xlabel: p.xlabel.clone(),
                    ylabel: p.ylabel.clone(),
                })
                .collect(),
        }
    }
    pub fn apply(&self, fig: &mut Figure) -> Result<()> {
        ensure!(
            (1..=4).contains(&self.columns) && self.columns <= self.panels.len().max(1),
            "invalid project columns"
        );
        let mut ids: Vec<_> = self.panels.iter().map(|p| p.id).collect();
        ids.sort_unstable();
        let mut expected: Vec<_> = fig.plots.iter().map(|p| p.id).collect();
        expected.sort_unstable();
        ensure!(
            ids == expected,
            "project panel identities do not match its recipe"
        );
        for text in [&self.title, &self.caption].into_iter().chain(
            self.panels
                .iter()
                .flat_map(|p| [&p.title, &p.caption, &p.xlabel, &p.ylabel]),
        ) {
            ensure!(text.len() <= 16_384, "project text exceeds 16384 bytes");
        }
        let mut plots: BTreeMap<_, _> = fig.plots.drain(..).map(|p| (p.id, p)).collect();
        fig.plots = self
            .panels
            .iter()
            .map(|p| {
                let mut plot = plots.remove(&p.id).expect("validated panel identity");
                plot.title = p.title.clone();
                plot.caption = p.caption.clone();
                plot.xlabel = p.xlabel.clone();
                plot.ylabel = p.ylabel.clone();
                plot
            })
            .collect();
        if let Some(subtitle) = &self.subtitle {
            fig.subtitle = subtitle.clone();
        }
        fig.title = self.title.clone();
        fig.caption = self.caption.clone();
        fig.columns = self.columns;
        Ok(())
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    pub presentation: Presentation,
    pub views: Vec<Bounds>,
    pub width: u32,
    pub height: u32,
    pub full_plot: bool,
}
impl Snapshot {
    pub fn validate(&self, fig: &Figure) -> Result<()> {
        ensure!(
            self.views.len() == fig.plots.len() && self.views.iter().all(|v| v.valid()),
            "invalid project views"
        );
        ensure!(
            (400..=8192).contains(&self.width)
                && (if fig.plots.first().is_some_and(|p| p.qc.margin_off) {
                    32
                } else {
                    300
                }..=8192)
                    .contains(&self.height)
                && self.width as u64 * self.height as u64 <= 32_000_000,
            "invalid project image size"
        );
        Ok(())
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Source {
    argument: usize,
    fingerprint: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Origin {
    arguments: Vec<String>,
    sources: Vec<Source>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    imports: Vec<ImportLayer>,
}
/// Source identity is the layer ID; series identity is (layer ID, original column index).
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Roles {
    pub columns: Vec<String>,
    pub x: Option<usize>,
    pub y: Vec<usize>,
    pub kind: String,
    pub missing: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ImportLayer {
    id: usize,
    origin: Box<Origin>,
    roles: Roles,
}
// FNV-1a is a deterministic change detector, not an authenticity/security signature.
pub fn fingerprint(path: &Path) -> Result<String> {
    let mut file =
        std::fs::File::open(path).with_context(|| format!("reading source {}", path.display()))?;
    let mut hash = 0xcbf29ce484222325u64;
    let mut size = 0u64;
    let mut buffer = [0; 65536];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        size += n as u64;
        for byte in &buffer[..n] {
            hash = (hash ^ *byte as u64).wrapping_mul(0x100000001b3);
        }
    }
    Ok(format!("fnv1a64:{size}:{hash:016x}"))
}
fn split_spec(spec: &str) -> (&str, &str) {
    let mut name = spec.strip_suffix('\'').unwrap_or(spec);
    while name.ends_with(']') || name.ends_with('}') {
        let open = if name.ends_with(']') { '[' } else { '{' };
        if let Some(at) = name.rfind(open) {
            name = &name[..at];
        } else {
            break;
        }
    }
    (name, &spec[name.len()..])
}
impl Origin {
    pub fn empty() -> Self {
        Self {
            arguments: vec![],
            sources: vec![],
            imports: vec![],
        }
    }
    pub fn add(&mut self, id: usize, origin: Origin, roles: Roles) {
        self.imports.push(ImportLayer {
            id,
            origin: Box::new(origin),
            roles,
        });
    }
    pub fn retain_imports(&mut self, ids: &std::collections::BTreeSet<usize>) {
        self.imports.retain(|layer| ids.contains(&layer.id));
    }
    fn verify_all(&self) -> Result<()> {
        self.verify()?;
        ensure!(self.imports.len() <= 64, "too many imported sources");
        let mut ids = std::collections::BTreeSet::new();
        for layer in &self.imports {
            ensure!(
                ids.insert(layer.id) && layer.origin.imports.is_empty(),
                "duplicate source identity or nested imports"
            );
            let role = &layer.roles;
            ensure!(
                !role.columns.is_empty()
                    && role.columns.len() <= 256
                    && !role.y.is_empty()
                    && role
                        .y
                        .iter()
                        .chain(role.x.iter())
                        .all(|&c| c < role.columns.len()),
                "invalid imported column roles"
            );
            layer.origin.verify()?;
        }
        Ok(())
    }
    fn rebuild(&self) -> Result<Figure> {
        self.verify_all()?;
        let mut fig = if self.arguments.is_empty() {
            Figure::empty()
        } else {
            model::from_options(&self.verify()?)?
        };
        let mut ids: std::collections::BTreeSet<_> = fig.plots.iter().map(|p| p.id).collect();
        for layer in &self.imports {
            ensure!(
                ids.insert(layer.id),
                "imported panel identity collides with original panel"
            );
            let mut imported = model::from_options(&layer.origin.verify()?)?;
            ensure!(
                imported.plots.len() == 1 && !imported.plots[0].qc.margin_off,
                "each imported recipe must produce one ordinary graph"
            );
            let mut plot = imported.plots.remove(0);
            plot.id = layer.id;
            fig.plots.push(plot);
            fig.synthetic |= imported.synthetic;
        }
        ensure!(fig.plots.len() <= 64, "project exceeds 64 graphs");
        ensure!(
            fig.plots.len() <= 1 || fig.plots.iter().all(|p| !p.qc.margin_off),
            "margin-free strips cannot be combined with other graphs"
        );
        Ok(fig)
    }
    fn rebase(&mut self, base: &Path, saving: bool) -> Result<()> {
        for source in &self.sources {
            let spec = self
                .arguments
                .get_mut(source.argument)
                .context("invalid source argument index")?;
            let (name, suffix) = split_spec(spec);
            let path = if saving {
                relative(Path::new(name), base)
            } else {
                base.join(name)
            };
            *spec = format!(
                "{}{suffix}",
                path.to_str().context("source path must be UTF-8")?
            );
        }
        for layer in &mut self.imports {
            ensure!(
                layer.origin.imports.is_empty(),
                "nested imports are unsupported"
            );
            layer.origin.rebase(base, saving)?;
        }
        Ok(())
    }
    pub fn capture(o: &cli::Options) -> Result<Self> {
        let mut origin = Self {
            arguments: crate::rsa::freeze(o)?,
            sources: vec![],
            imports: vec![],
        };
        // Unit-test/demo callers may construct Options directly.
        if origin.arguments.is_empty() && o.files.is_empty() {
            origin.arguments = vec!["-demo".into()];
        }
        for &argument in &o.source_args {
            let spec = &origin.arguments[argument];
            ensure!(
                spec != "-" && spec != "-stdin" && !split_spec(spec).0.eq("-"),
                "Save Project needs reusable input: save stdin to a file and reopen it first"
            );
            if spec.starts_with("1D:") {
                continue;
            }
            let (name, suffix) = split_spec(spec);
            let source = Path::new(name);
            let parent = source
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            // Preserve symlink basenames: automatic series labels use them too.
            let path = std::fs::canonicalize(parent)
                .with_context(|| format!("locating source {name}"))?
                .join(source.file_name().context("source needs a filename")?);
            let hash = fingerprint(&path)?;
            origin.arguments[argument] = format!(
                "{}{suffix}",
                path.to_str().context("source path must be UTF-8")?
            );
            origin.sources.push(Source {
                argument,
                fingerprint: hash,
            });
        }
        Ok(origin)
    }
    pub fn verify(&self) -> Result<cli::Options> {
        if self.arguments.is_empty() {
            ensure!(
                self.sources.is_empty(),
                "empty recipe has source references"
            );
            return Ok(cli::Options::default());
        }
        let o = cli::parse(&self.arguments).context("invalid project recipe")?;
        let expected: Vec<_> = o
            .source_args
            .iter()
            .copied()
            .filter(|&i| !self.arguments[i].starts_with("1D:"))
            .collect();
        ensure!(
            expected == self.sources.iter().map(|s| s.argument).collect::<Vec<_>>(),
            "project source manifest does not match recipe"
        );
        for s in &self.sources {
            let name = split_spec(&self.arguments[s.argument]).0;
            ensure!(
                name != "-" && name != "-stdin",
                "projects cannot replay stdin"
            );
            ensure!(
                fingerprint(Path::new(name))? == s.fingerprint,
                "source changed: {name}. Restore the original data or start a new figure from the changed source"
            );
        }
        Ok(o)
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Project {
    format: String,
    version: u32,
    origin: Origin,
    state: Snapshot,
}
pub struct Loaded {
    pub figure: Figure,
    pub origin: Origin,
    pub state: Snapshot,
}
fn relative(path: &Path, base: &Path) -> PathBuf {
    let p: Vec<_> = path.components().collect();
    let b: Vec<_> = base.components().collect();
    let n = p.iter().zip(&b).take_while(|(a, b)| a == b).count();
    if n == 0 {
        return path.to_path_buf();
    }
    let mut out = PathBuf::new();
    for _ in n..b.len() {
        out.push("..");
    }
    for c in &p[n..] {
        out.push(c);
    }
    out
}
pub fn save(path: &Path, origin: &Origin, state: &Snapshot, overwrite: bool) -> Result<()> {
    let mut recipe = origin.clone();
    recipe.retain_imports(&state.presentation.panels.iter().map(|p| p.id).collect());
    recipe.verify_all()?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let base = std::fs::canonicalize(parent)?;
    recipe.rebase(&base, true)?;
    let project = Project {
        format: "1dPetePlot".into(),
        version: 2,
        origin: recipe,
        state: state.clone(),
    };
    let mut temp = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut temp, &project)?;
    temp.write_all(b"\n")?;
    temp.flush()?;
    if overwrite {
        temp.persist(path)?;
    } else {
        temp.persist_noclobber(path)?;
    }
    Ok(())
}
pub fn load(path: &Path) -> Result<Loaded> {
    let file =
        std::fs::File::open(path).with_context(|| format!("opening project {}", path.display()))?;
    ensure!(
        file.metadata()?.len() <= 16_000_000,
        "project exceeds 16 MB recipe limit"
    );
    let mut p: Project = serde_json::from_reader(file).context("reading project JSON")?;
    ensure!(
        p.format == "1dPetePlot" && (1..=2).contains(&p.version),
        "unsupported project format/version (expected 1dPetePlot version 1 or 2)"
    );
    let base = std::fs::canonicalize(path)?.parent().unwrap().to_path_buf();
    p.origin.rebase(&base, false)?;
    let mut figure = p.origin.rebuild().context("rebuilding project data")?;
    p.state.presentation.apply(&mut figure)?;
    p.state.validate(&figure)?;
    Ok(Loaded {
        figure,
        origin: p.origin,
        state: p.state,
    })
}

#[cfg_attr(not(feature = "gui"), allow(dead_code))]
#[derive(Default)]
pub struct History {
    undo: Vec<Presentation>,
    redo: Vec<Presentation>,
}
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
impl History {
    pub fn record(&mut self, before: Presentation, after: &Presentation) {
        if &before == after {
            return;
        }
        self.undo.push(before);
        self.redo.clear();
        if self.undo.len() > 100 {
            self.undo.remove(0);
        }
    }
    pub fn undo(&mut self, current: Presentation) -> Option<Presentation> {
        let previous = self.undo.pop()?;
        self.redo.push(current);
        Some(previous)
    }
    pub fn redo(&mut self, current: Presentation) -> Option<Presentation> {
        let next = self.redo.pop()?;
        self.undo.push(current);
        Some(next)
    }
    pub fn panel_ids(&self) -> std::collections::BTreeSet<usize> {
        self.undo
            .iter()
            .chain(&self.redo)
            .flat_map(|p| p.panels.iter().map(|p| p.id))
            .collect()
    }
    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn example(args: &[String]) -> (Origin, Figure, Snapshot) {
        let o = cli::parse(args).unwrap();
        let fig = model::from_options(&o).unwrap();
        let state = Snapshot {
            presentation: Presentation::capture(&fig),
            views: fig.plots.iter().map(|p| p.bounds).collect(),
            width: 1200,
            height: 1000,
            full_plot: false,
        };
        (Origin::capture(&o).unwrap(), fig, state)
    }
    #[test]
    fn project_roundtrip_relocates_sources_preserves_missing_values_and_pixels() {
        let temp = tempfile::tempdir().unwrap();
        let old = temp.path().join("old");
        std::fs::create_dir(&old).unwrap();
        let input = old.join("values.tsv");
        std::fs::write(&input, "a\tb\n1\t4\nNaN\t3\n2\t6\n").unwrap();
        #[cfg(unix)]
        let input = {
            let alias = old.join("alias.tsv");
            std::os::unix::fs::symlink("values.tsv", &alias).unwrap();
            alias
        };
        let (origin, mut fig, mut state) = example(&[
            "-infiles".into(),
            format!("{}[a,b]", input.display()),
            "-missing".into(),
            "gap".into(),
            "-legend_on".into(),
        ]);
        state.presentation.title = "Effects α < β".into();
        state.presentation.caption = "Figure caption with scientific r^{2} and two lines.\nUnchanged sources, reproducible edits.".into();
        state.presentation.panels[0].caption = "Missing observations remain missing.".into();
        state.presentation.panels.reverse();
        state.views.reverse();
        state.views[0] = state.views[0].pan(0.2, 0.1);
        state.presentation.apply(&mut fig).unwrap();
        let path = old.join("figure.peteplot.json");
        save(&path, &origin, &state, false).unwrap();
        assert!(save(&path, &origin, &state, false).is_err());
        let json = std::fs::read_to_string(&path).unwrap();
        assert!(!json.contains(old.to_str().unwrap()));
        assert!(json.contains(&format!(
            "{}[a,b]",
            input.file_name().unwrap().to_string_lossy()
        )));
        let new = temp.path().join("moved");
        std::fs::rename(&old, &new).unwrap();
        let reopened = load(&new.join("figure.peteplot.json")).unwrap();
        assert_eq!(reopened.state, state);
        let model::Content::Lines(ss) = &reopened.figure.plots[1].content else {
            panic!()
        };
        assert!(ss[0].points[1][1].is_nan());
        let renderer = crate::render::Renderer::new();
        let render = |f: &Figure| {
            renderer
                .raster(
                    &renderer
                        .scene(f, &state.views, state.width, state.height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(render(&fig).data(), render(&reopened.figure).data());
        // A changed or missing source is never silently accepted.
        std::fs::write(new.join("values.tsv"), "a\tb\n1\t4\n8\t3\n2\t6\n").unwrap();
        assert!(
            format!(
                "{:#}",
                load(&new.join("figure.peteplot.json")).err().unwrap()
            )
            .contains("source changed")
        );
        std::fs::remove_file(new.join("values.tsv")).unwrap();
        assert!(load(&new.join("figure.peteplot.json")).is_err());
    }
    #[test]
    fn validates_version_identity_views_manifest_and_detects_changes_before_save() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("figure.json");
        let (origin, _, state) = example(&["-demo".into()]);
        save(&path, &origin, &state, false).unwrap();
        let original: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        for (field, value) in [
            ("/version", serde_json::json!(99)),
            ("/state/presentation/columns", serde_json::json!(0)),
            ("/state/presentation/panels/0/id", serde_json::json!(99)),
            ("/state/views/0/x", serde_json::json!([0, 0])),
            ("/state/width", serde_json::json!(90000)),
        ] {
            let mut bad = original.clone();
            *bad.pointer_mut(field).unwrap() = value;
            std::fs::write(&path, serde_json::to_vec(&bad).unwrap()).unwrap();
            assert!(load(&path).is_err(), "{field}");
        }
        let file = temp.path().join("source.1D");
        std::fs::write(&file, "1\n2\n3\n").unwrap();
        let (origin, _, state) = example(&[file.display().to_string()]);
        std::fs::write(file, "4\n5\n6\n").unwrap();
        assert!(save(&path, &origin, &state, true).is_err());
    }
    #[test]
    fn tracks_only_source_arguments_and_all_auxiliary_inputs() {
        let temp = tempfile::tempdir().unwrap();
        let file = temp.path().join("source.1D");
        std::fs::write(&file, "1\n2\n3\n").unwrap();
        let name = file.display().to_string();
        let o = cli::parse(
            &[
                "-title",
                &name,
                "-xfile",
                &name,
                "-yfiles_pm",
                &name,
                "-censor_files",
                &name,
                "-shade-xmat",
                &name,
                "-infiles",
                &name,
            ]
            .map(str::to_owned),
        )
        .unwrap();
        assert_eq!(o.source_args, vec![3, 5, 7, 9, 11]);
        let origin = Origin::capture(&o).unwrap();
        assert_eq!(origin.arguments[1], name);
        assert_eq!(origin.sources.len(), 5);
        for arg in ["-stdin", "-", "-'", "-[0]"] {
            let args = if arg == "-stdin" {
                vec![arg.into()]
            } else {
                vec!["--".into(), arg.into()]
            };
            let o = cli::parse(&args).unwrap();
            assert!(Origin::capture(&o).is_err());
        }
    }
    #[test]
    fn inline_recipe_and_history_branching() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("inline.json");
        let (origin, _, state) =
            example(&["1D: 0 NaN 2 3".into(), "-missing".into(), "gap".into()]);
        save(&path, &origin, &state, false).unwrap();
        assert!(load(&path).is_ok());
        let mut history = History::default();
        let before = state.presentation;
        let mut after = before.clone();
        after.caption = "New".into();
        history.record(before.clone(), &after);
        assert_eq!(history.undo(after.clone()), Some(before.clone()));
        assert_eq!(history.redo(before.clone()), Some(after.clone()));
        history.undo(after);
        let mut branch = before.clone();
        branch.title = "New branch".into();
        history.record(before, &branch);
        assert!(!history.can_redo());
    }
}
