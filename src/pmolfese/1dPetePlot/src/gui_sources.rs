//! Explicit source recovery/reload. Background candidates cannot mutate the open document.
use crate::{
    document::{self, Loaded, Origin, Snapshot},
    model::Figure,
    render::Renderer,
};
use anyhow::{Context, Result, ensure};
use eframe::egui;
use std::{
    path::{Path, PathBuf},
    sync::mpsc,
};

pub struct Preview {
    pub project: Loaded,
    pub rgba: Vec<u8>,
    pub size: [usize; 2],
    pub report: Vec<String>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Choice {
    pub original: String,
    pub replacement: String,
    pub status: String,
}
pub fn prepare(
    origin: &Origin,
    state: &Snapshot,
    before: Option<&Figure>,
    choices: &[Choice],
    keep_views: bool,
) -> Result<Preview> {
    let replacements = choices
        .iter()
        .map(|c| (c.original.clone(), c.replacement.clone()))
        .collect();
    for c in choices {
        let meta = std::fs::metadata(&c.replacement)
            .with_context(|| format!("opening {}", c.replacement))?;
        ensure!(
            meta.is_file() && meta.len() <= 64 * 1024 * 1024,
            "source reload supports regular files up to 64 MiB each"
        );
    }
    let next = origin.refreshed(&replacements)?;
    // A same-sized reordered table is not the same column mapping.
    for (o, roles) in next.role_recipes() {
        let file = o
            .rsa
            .file
            .as_ref()
            .or(o.rdm.files.first())
            .or(o.files.first())
            .context("import recipe has no table")?;
        let table = crate::import::inspect(Path::new(document::split_spec(file).0), || false)?;
        ensure!(
            roles.columns
                == table
                    .columns
                    .iter()
                    .map(|c| c.name.clone())
                    .collect::<Vec<_>>(),
            "column names/order changed in {}; reimport to assign roles explicitly",
            file
        );
    }
    let mut figure = next.rebuild()?;
    let mut report = vec![];
    if let Some(old) = before {
        ensure!(
            old.plots.len() == figure.plots.len(),
            "reload changed the number of graphs; reimport this result to assign panels explicitly"
        );
        for a in &old.plots {
            let b = figure
                .plots
                .iter()
                .find(|p| p.id == a.id)
                .context("graph identity changed")?;
            ensure!(
                std::mem::discriminant(&a.content) == std::mem::discriminant(&b.content)
                    && a.style_series().len() == b.style_series().len(),
                "graph {} changed plot/series roles; reimport to assign roles",
                a.id + 1
            );
            report.push(format!(
                "{}: {} → {} · range {:?} → {:?}",
                a.title,
                observations(a),
                observations(b),
                a.bounds.y,
                b.bounds.y
            ));
            let old_labels: Vec<_> = a.style_series().iter().map(|s| s.0.to_string()).collect();
            let new_labels: Vec<_> = b.style_series().iter().map(|s| s.0.to_string()).collect();
            if old_labels != new_labels {
                report.push(format!(
                    "Series labels: {} → {}. Review retained legend overrides.",
                    old_labels.join(", "),
                    new_labels.join(", ")
                ));
            }
        }
    }
    let mut state = state.clone();
    state.presentation.apply(&mut figure)?;
    if !keep_views {
        state.views = figure.plots.iter().map(|p| p.bounds).collect();
    }
    state.height = state
        .height
        .max(crate::rsa::minimum_height(&figure, state.width));
    state.validate(&figure)?;
    let renderer = Renderer::new();
    renderer.scene(&figure, &state.views, state.width, state.height)?;
    let width = 900;
    let height = (state.height as f64 * width as f64 / state.width as f64)
        .round()
        .max(300.) as u32;
    let scene = renderer.scene(&figure, &state.views, width, height)?;
    let pix = renderer.raster(&scene)?;
    next.verify_all()?;
    Ok(Preview {
        project: Loaded {
            figure,
            origin: next,
            state,
        },
        rgba: pix.data().to_vec(),
        size: [width as usize, height as usize],
        report,
    })
}
fn observations(p: &crate::model::Plot) -> usize {
    use crate::model::Content;
    match &p.content {
        Content::Lines(s) | Content::Scatter(s) => s.iter().map(|s| s.points.len()).sum(),
        Content::Histogram { series, .. } => series.iter().map(|s| s.points.len()).sum(),
        Content::Intervals { rows, .. } => rows.len(),
        Content::Heat { matrix, .. } => matrix.iter().map(Vec::len).sum(),
    }
}
enum Reply {
    Status(Vec<String>),
    Preview(Box<Result<Preview, String>>, Vec<Choice>, bool),
}
pub struct Editor {
    origin: Origin,
    state: Snapshot,
    before: Option<Figure>,
    pub recovery: Option<PathBuf>,
    choices: Vec<Choice>,
    keep_views: bool,
    receiver: Option<mpsc::Receiver<Reply>>,
    candidate: Option<Preview>,
    texture: Option<egui::TextureHandle>,
    error: Option<String>,
    #[cfg(test)]
    pub apply_rect: Option<egui::Rect>,
}
impl Editor {
    pub fn new(
        origin: Origin,
        state: Snapshot,
        before: Option<Figure>,
        recovery: Option<PathBuf>,
        ctx: &egui::Context,
    ) -> Self {
        let paths = origin.source_paths();
        let choices = paths
            .iter()
            .map(|p| Choice {
                original: p.clone(),
                replacement: p.clone(),
                status: "Checking…".into(),
            })
            .collect();
        let (tx, rx) = mpsc::channel();
        let source = origin.clone();
        let context = ctx.clone();
        std::thread::spawn(move || {
            let statuses = paths.iter().map(|p| source.source_status(p)).collect();
            let _ = tx.send(Reply::Status(statuses));
            context.request_repaint();
        });
        Self {
            origin,
            state,
            before,
            recovery,
            choices,
            keep_views: false,
            receiver: Some(rx),
            candidate: None,
            texture: None,
            error: None,
            #[cfg(test)]
            apply_rect: None,
        }
    }
    pub fn reject(&mut self, message: String) {
        self.error = Some(message);
        self.candidate = None;
        self.texture = None;
    }
    fn preview(&mut self, ctx: &egui::Context) {
        self.candidate = None;
        self.texture = None;
        self.error = None;
        let origin = self.origin.clone();
        let state = self.state.clone();
        let before = self.before.clone();
        let choices = self.choices.clone();
        let keep = self.keep_views;
        let ctx = ctx.clone();
        let (tx, rx) = mpsc::channel();
        self.receiver = Some(rx);
        std::thread::spawn(move || {
            let answer = prepare(&origin, &state, before.as_ref(), &choices, keep)
                .map_err(|e| format!("{e:#}"));
            let _ = tx.send(Reply::Preview(Box::new(answer), choices, keep));
            ctx.request_repaint();
        });
    }
    /// Returns (cancelled, ready-to-apply candidate); caller rechecks fingerprints before mutation.
    pub fn show(&mut self, ui: &mut egui::Ui) -> (bool, Option<Loaded>) {
        let ctx = ui.ctx().clone();
        if let Some(reply) = self.receiver.as_ref().and_then(|r| r.try_recv().ok()) {
            self.receiver = None;
            match reply {
                Reply::Status(statuses) => {
                    for (c, status) in self.choices.iter_mut().zip(statuses) {
                        c.status = status;
                    }
                }
                Reply::Preview(result, choices, keep)
                    if choices == self.choices && keep == self.keep_views =>
                {
                    match *result {
                        Ok(p) => {
                            self.texture = Some(ctx.load_texture(
                                "source-preview",
                                egui::ColorImage::from_rgba_premultiplied(p.size, &p.rgba),
                                egui::TextureOptions::LINEAR,
                            ));
                            self.candidate = Some(p);
                        }
                        Err(e) => self.error = Some(e),
                    }
                }
                _ => {}
            }
        }
        let mut cancel = false;
        let mut apply = false;
        let mut preview = false;
        let mut changed = false;
        egui::CentralPanel::default().show(ui,|ui|{
            let height=(ui.available_height()-48.).max(100.);
            egui::ScrollArea::vertical().max_height(height).auto_shrink([false,false]).show(ui,|ui|{
                ui.heading(if self.recovery.is_some(){"Recover project sources"}else{"Reload / relink sources"});
                ui.label("Choose replacements, preview the data changes, then Apply. Titles, captions and compatible styling are retained; review their scientific meaning after changing data.");
                for c in &mut self.choices {
                    ui.separator();ui.label(&c.original);ui.small(&c.status);
                    ui.horizontal(|ui|{changed|=ui.text_edit_singleline(&mut c.replacement).changed();if ui.button("Choose…").clicked() && let Some(path)=rfd::FileDialog::new().set_title("Relink source").pick_file(){c.replacement=path.to_string_lossy().into();changed=true;}});
                }
                if self.choices.is_empty(){ui.label("This figure has no external source files.");}
                changed|=ui.checkbox(&mut self.keep_views,"Keep current axis views (otherwise fit refreshed data)").changed();
                if let Some(p)=&self.candidate {for line in &p.report {ui.small(line);}}
                if let Some(t)=&self.texture {let size=t.size_vec2();let width=ui.available_width().min(size.x);ui.add(egui::Image::new((t.id(),egui::vec2(width,width*size.y/size.x))));}
                if let Some(e)=&self.error {ui.colored_label(egui::Color32::DARK_RED,e);}
                if self.receiver.is_some(){ui.spinner();ctx.request_repaint_after(std::time::Duration::from_millis(50));}
            });
            ui.separator();ui.horizontal(|ui|{
                cancel=ui.button("Cancel").clicked();
                preview=ui.add_enabled(self.receiver.is_none() && !self.choices.is_empty(),egui::Button::new("Preview changes")).clicked();
                let response=ui.add_enabled(self.candidate.is_some() && !changed && self.receiver.is_none(),egui::Button::new("Apply changes"));apply=response.clicked();
                #[cfg(test)] {self.apply_rect=Some(response.rect);}
            });
        });
        if changed {
            self.candidate = None;
            self.texture = None;
        }
        if preview {
            self.preview(&ctx);
        }
        if apply && let Some(candidate) = self.candidate.take() {
            match candidate.project.origin.verify_all() {
                Ok(()) => return (false, Some(candidate.project)),
                Err(e) => {
                    self.error = Some(format!(
                        "Source changed after preview; preview again: {e:#}"
                    ));
                    self.texture = None;
                }
            }
        }
        (cancel, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn setup(dir: &Path) -> (Origin, Snapshot, Figure, String) {
        let p = dir.join("values.tsv");
        std::fs::write(&p, "time\teffect\n0\t1\n1\t2\n2\t3\n").unwrap();
        let table = crate::import::inspect(&p, || false).unwrap();
        let mut mapping = crate::import::Mapping::for_kind(&table, crate::import::Kind::Line);
        mapping.x = Some(0);
        mapping.y = vec![1];
        let prepared = crate::import::prepare(&table, &mapping, || false).unwrap();
        let mut origin = Origin::empty();
        origin.add(0, prepared.origin, prepared.roles);
        let f = prepared.figure;
        let state = Snapshot {
            presentation: document::Presentation::capture(&f),
            views: vec![f.plots[0].bounds],
            width: 1000,
            height: 700,
            full_plot: false,
        };
        {
            let file = origin.source_paths()[0].clone();
            (origin, state, f, file)
        }
    }
    #[test]
    fn refresh_relink_and_recovery_validate_roles_and_never_change_original_recipe() {
        let d = tempfile::tempdir().unwrap();
        let (origin, state, f, file) = setup(d.path());
        let project = d.path().join("figure.json");
        document::save(&project, &origin, &state, false).unwrap();
        std::fs::write(&file, "time\teffect\n0\t4\n1\t5\n2\t6\n").unwrap();
        assert_eq!(origin.source_status(&file), "Changed");
        assert!(document::load(&project).is_err());
        let choices = vec![Choice {
            original: file.clone(),
            replacement: file.clone(),
            status: "Changed".into(),
        }];
        let p = prepare(&origin, &state, Some(&f), &choices, false).unwrap();
        assert!(origin.verify_all().is_err());
        p.project.origin.verify_all().unwrap();
        let crate::model::Content::Lines(ss) = &p.project.figure.plots[0].content else {
            panic!()
        };
        assert_eq!(ss[0].points[0], [0., 4.]);
        let relocated = d.path().join("relocated.tsv");
        std::fs::rename(&file, &relocated).unwrap();
        assert!(origin.source_status(&file).starts_with("Missing"));
        let choices = vec![Choice {
            original: file.clone(),
            replacement: relocated.to_string_lossy().into(),
            status: "Missing".into(),
        }];
        let (raw, state) = document::read_project(&project).unwrap();
        let p = prepare(&raw, &state, None, &choices, false).unwrap();
        document::save(&project, &p.project.origin, &p.project.state, true).unwrap();
        document::load(&project).unwrap();
        std::fs::write(&relocated, "effect\ttime\n4\t0\n5\t1\n").unwrap();
        assert!(
            prepare(&origin, &state, Some(&f), &choices, false)
                .err()
                .unwrap()
                .to_string()
                .contains("column names/order")
        );
    }
    #[test]
    fn source_review_footer_is_visible_and_stale_candidate_is_refused() {
        let d = tempfile::tempdir().unwrap();
        let (origin, state, f, file) = setup(d.path());
        let ctx = egui::Context::default();
        let mut editor = Editor::new(origin, state, Some(f), None, &ctx);
        let frame = |editor: &mut Editor| {
            let mut result = None;
            ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(egui::Rect::from_min_size(
                        egui::Pos2::ZERO,
                        egui::vec2(850., 500.),
                    )),
                    ..Default::default()
                },
                |ui| result = Some(editor.show(ui)),
            )
            .drop_without_applying_deltas();
            result.unwrap()
        };
        let wait = |editor: &mut Editor| {
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
            while editor.receiver.is_some() {
                frame(editor);
                assert!(std::time::Instant::now() < deadline);
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
        };
        wait(&mut editor);
        editor.preview(&ctx);
        wait(&mut editor);
        assert!(editor.candidate.is_some(), "{:?}", editor.error);
        frame(&mut editor);
        assert!(ctx.content_rect().contains_rect(editor.apply_rect.unwrap()));
        std::fs::write(&file, "time\teffect\n0\t99\n").unwrap();
        let pos = editor.apply_rect.unwrap().center();
        let mut result = (false, None);
        for pressed in [true, false] {
            ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(ctx.content_rect()),
                    events: vec![
                        egui::Event::PointerMoved(pos),
                        egui::Event::PointerButton {
                            pos,
                            button: egui::PointerButton::Primary,
                            pressed,
                            modifiers: egui::Modifiers::NONE,
                        },
                    ],
                    ..Default::default()
                },
                |ui| result = editor.show(ui),
            )
            .drop_without_applying_deltas();
        }
        assert!(result.1.is_none());
        assert!(
            editor
                .error
                .as_ref()
                .unwrap()
                .contains("Source changed after preview")
        );
    }
}
