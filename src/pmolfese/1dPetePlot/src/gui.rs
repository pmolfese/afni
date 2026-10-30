#[path = "gui_sources.rs"]
mod sources;
#[path = "gui_style.rs"]
mod styling;
use crate::{
    cli::Options,
    document,
    model::{Bounds, Figure, Plot},
    render::{Area, Renderer, TextRegion, TextTarget},
};
use anyhow::{Result, anyhow};
use eframe::egui::{self, Color32, RichText};
use std::path::PathBuf;

pub fn run(
    fig: Figure,
    o: Options,
    renderer: Renderer,
    loaded: Option<document::Loaded>,
    project_path: Option<PathBuf>,
    origin: Option<Result<document::Origin>>,
) -> Result<()> {
    #[cfg(target_os = "linux")]
    if std::env::var_os("WAYLAND_DISPLAY").is_none() {
        return Err(anyhow!(
            "native viewing requires Wayland on Linux; use -prefix plot.png for headless export"
        ));
    }
    let native = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1320., 900.])
            .with_min_inner_size([850., 600.])
            .with_decorations(true)
            .with_title("1dPetePlot")
            .with_icon(window_icon()?),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    eframe::run_native(
        "1dPetePlot",
        native,
        Box::new(move |cc| {
            cc.egui_ctx.set_visuals(egui::Visuals::light());
            let mut style = (*cc.egui_ctx.style_of(egui::Theme::Light)).clone();
            style.spacing.item_spacing = egui::vec2(10., 10.);
            style.visuals.selection.bg_fill = Color32::from_rgb(20, 125, 146);
            style.visuals.selection.stroke.color = Color32::WHITE;
            cc.egui_ctx.set_style_of(egui::Theme::Light, style);
            let dimensions = (o.width, o.height);
            let mut app = App::new(fig, o, renderer);
            if let Some(origin) = origin {
                match origin {
                    Ok(origin) => {
                        app.origin = Some(origin);
                        app.origin_error.clear();
                    }
                    Err(error) => {
                        app.origin = None;
                        app.origin_error = format!("Cannot save project: {error:#}");
                    }
                }
            }
            if let Some(project) = loaded {
                app.install_project(project, project_path.expect("recipe path"));
                app.export_width = dimensions.0;
                app.export_height = dimensions.1;
            }
            #[cfg(target_os = "macos")]
            let app = {
                let mut app = app;
                app.menu = Some(crate::menus::NativeMenu::new(&cc.egui_ctx)?);
                app
            };
            Ok(Box::new(app))
        }),
    )
    .map_err(|e| anyhow!("opening native window: {e}"))
}
fn is_fullscreen(ctx: &egui::Context) -> bool {
    ctx.input(|i| i.viewport().fullscreen.unwrap_or(false))
}
fn exit_fullscreen_on_escape(ctx: &egui::Context) {
    if is_fullscreen(ctx)
        && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape))
    {
        ctx.send_viewport_cmd(egui::ViewportCommand::Fullscreen(false));
    }
}
fn window_icon() -> Result<egui::IconData> {
    let tree =
        resvg::usvg::Tree::from_str(include_str!("../assets/peteplot.svg"), &Default::default())?;
    let mut pix =
        resvg::tiny_skia::Pixmap::new(256, 256).ok_or_else(|| anyhow!("allocating icon"))?;
    resvg::render(&tree, Default::default(), &mut pix.as_mut());
    Ok(egui::IconData {
        rgba: pix
            .pixels()
            .iter()
            .flat_map(|p| {
                let c = p.demultiply();
                [c.red(), c.green(), c.blue(), c.alpha()]
            })
            .collect(),
        width: 256,
        height: 256,
    })
}
#[derive(Clone)]
struct SaveRequest {
    fig: Figure,
    views: Vec<Bounds>,
    width: u32,
    height: u32,
    path: PathBuf,
}
struct TextEditor {
    target: TextTarget,
    value: String,
    first: bool,
    position: Option<egui::Pos2>,
}
#[derive(Clone)]
struct DataVersion {
    width: u32,
    height: u32,
    full_plot: bool,
    figure: Figure,
    origin: Option<document::Origin>,
    views: Vec<Bounds>,
    dormant: std::collections::BTreeMap<usize, Plot>,
    dormant_views: std::collections::BTreeMap<usize, Bounds>,
}
struct App {
    sources: Option<sources::Editor>,
    data_revision: u64,
    saved_data_revision: u64,
    data_versions: std::collections::BTreeMap<u64, DataVersion>,
    import: crate::import_ui::ImportUi,
    dormant: std::collections::BTreeMap<usize, Plot>,
    dormant_views: std::collections::BTreeMap<usize, Bounds>,
    next_id: usize,
    origin: Option<document::Origin>,
    origin_error: String,
    project_path: Option<PathBuf>,
    saved: document::Snapshot,
    history: document::History,
    editor: Option<TextEditor>,
    styling: Option<styling::Editor>,
    duplicate_pending: Option<usize>,
    closing: bool,
    window_title: String,
    text_regions: Vec<TextRegion>,
    fig: Figure,
    views: Vec<Bounds>,
    focus: Option<usize>,
    renderer: Renderer,
    texture: Option<egui::TextureHandle>,
    thumbnails: Vec<egui::TextureHandle>,
    thumbnail_titles: Vec<String>,
    areas: Vec<Area>,
    logical_height: f64,
    dirty: bool,
    last_size: [u32; 2],
    drag_panel: Option<usize>,
    zoom_start: Option<egui::Pos2>,
    context_panel: Option<usize>,
    edit_panel: usize,
    clipboard: Option<arboard::Clipboard>,
    #[cfg(target_os = "macos")]
    menu: Option<crate::menus::NativeMenu>,
    save_path: String,
    export_width: u32,
    export_height: u32,
    full_plot: bool,
    pending: Option<SaveRequest>,
    status: String,
    hover: String,
    error: bool,
}
impl App {
    fn new(fig: Figure, o: Options, renderer: Renderer) -> Self {
        let views: Vec<_> = fig.plots.iter().map(|p| p.bounds).collect();
        let saved = document::Snapshot {
            presentation: document::Presentation::capture(&fig),
            views: views.clone(),
            width: o.width,
            height: o.height,
            full_plot: false,
        };
        let next_id = fig.plots.iter().map(|p| p.id + 1).max().unwrap_or(0);
        let source = if fig.plots.is_empty() {
            Ok(document::Origin::empty())
        } else {
            document::Origin::capture(&o)
        };
        let (origin, origin_error) = match source {
            Ok(o) => (Some(o), String::new()),
            Err(e) => (None, format!("Cannot save project: {e:#}")),
        };
        App {
            sources: None,
            data_revision: 0,
            saved_data_revision: 0,
            data_versions: Default::default(),
            import: Default::default(),
            dormant: Default::default(),
            dormant_views: Default::default(),
            next_id,
            origin,
            origin_error,
            saved,
            project_path: None,
            history: Default::default(),
            editor: None,
            styling: None,
            duplicate_pending: None,
            closing: false,
            window_title: String::new(),
            text_regions: vec![],
            fig,
            views,
            focus: None,
            renderer,
            texture: None,
            thumbnails: vec![],
            thumbnail_titles: vec![],
            areas: vec![],
            logical_height: 0.,
            dirty: true,
            last_size: [0, 0],
            drag_panel: None,
            zoom_start: None,
            context_panel: None,
            edit_panel: 0,
            clipboard: None,
            #[cfg(target_os = "macos")]
            menu: None,
            save_path: o.output.unwrap_or_else(|| "1dPetePlot.png".into()),
            export_width: o.width,
            export_height: o.height,
            full_plot: false,
            pending: None,
            status: "Ready · drag to pan · Shift-drag to zoom · double-click to reset".into(),
            hover: String::new(),
            error: false,
        }
    }
    fn snapshot(&self) -> document::Snapshot {
        document::Snapshot {
            presentation: document::Presentation::capture(&self.fig),
            views: self.views.clone(),
            width: self.export_width,
            height: self.export_height,
            full_plot: self.full_plot,
        }
    }
    fn modified(&self) -> bool {
        self.editor.is_some()
            || self.styling.is_some()
            || self.snapshot() != self.saved
            || self.data_revision != self.saved_data_revision
    }
    fn edited(&mut self, before: document::Presentation) {
        self.history
            .record(before, &document::Presentation::capture(&self.fig));
        let mut ids = self.history.panel_ids();
        ids.extend(self.fig.plots.iter().map(|p| p.id));
        self.dormant.retain(|id, _| ids.contains(id));
        self.dormant_views.retain(|id, _| ids.contains(id));
        if let Some(origin) = &mut self.origin {
            origin.retain_imports(&ids);
        }
        let revisions = self.history.revisions();
        self.data_versions.retain(|r, _| revisions.contains(r));
        self.thumbnails.clear();
        self.dirty = true;
    }
    fn restore_presentation(&mut self, presentation: document::Presentation) {
        self.dormant_views.extend(
            self.fig
                .plots
                .iter()
                .zip(&self.views)
                .map(|(p, v)| (p.id, *v)),
        );
        let focus = self.focus.and_then(|i| self.fig.plots.get(i)).map(|p| p.id);
        let edit = self.fig.plots.get(self.edit_panel).map(|p| p.id);
        self.dormant
            .extend(self.fig.plots.drain(..).map(|p| (p.id, p)));
        self.fig.plots = presentation
            .panels
            .iter()
            .map(|p| {
                self.dormant
                    .remove(&p.id)
                    .expect("history owns its panel data")
            })
            .collect();
        presentation
            .apply(&mut self.fig)
            .expect("validated presentation history");
        self.views = self
            .fig
            .plots
            .iter()
            .map(|p| self.dormant_views.remove(&p.id).unwrap_or(p.bounds))
            .collect();
        self.focus = focus.and_then(|id| self.fig.plots.iter().position(|p| p.id == id));
        self.edit_panel = edit
            .and_then(|id| self.fig.plots.iter().position(|p| p.id == id))
            .unwrap_or(0);
        self.texture = None;
        self.areas.clear();
        self.text_regions.clear();
        self.drag_panel = None;
        self.zoom_start = None;
        self.context_panel = None;
        self.thumbnails.clear();
        self.dirty = true;
    }
    fn undo(&mut self, redo: bool) {
        if self.editor.is_some() || self.styling.is_some() {
            return;
        }
        let current = document::Presentation::capture(&self.fig);
        let next = if redo {
            self.history.redo(current)
        } else {
            self.history.undo(current)
        };
        if let Some(p) = next {
            if self.data_revision != self.history.revision() {
                self.cache_data();
                let version = self.data_versions[&self.history.revision()].clone();
                self.export_width = version.width;
                self.export_height = version.height;
                self.full_plot = version.full_plot;
                self.fig = version.figure;
                self.origin = version.origin;
                self.views = version.views;
                self.dormant = version.dormant;
                self.dormant_views = version.dormant_views;
                self.data_revision = self.history.revision();
            }
            self.restore_presentation(p);
            self.status = if redo { "Redone" } else { "Undone" }.into();
            self.error = false;
        }
    }
    fn cache_data(&mut self) {
        self.data_versions.insert(
            self.data_revision,
            DataVersion {
                width: self.export_width,
                height: self.export_height,
                full_plot: self.full_plot,
                figure: self.fig.clone(),
                origin: self.origin.clone(),
                views: self.views.clone(),
                dormant: self.dormant.clone(),
                dormant_views: self.dormant_views.clone(),
            },
        );
    }
    fn apply_data(&mut self, project: document::Loaded) -> Result<()> {
        project.origin.verify_all()?;
        project.state.validate(&project.figure)?;
        self.renderer.scene(
            &project.figure,
            &project.state.views,
            project.state.width,
            project.state.height,
        )?;
        self.cache_data();
        self.history
            .record_data(document::Presentation::capture(&self.fig));
        self.data_revision = self.history.revision();
        self.fig = project.figure;
        self.origin = Some(project.origin);
        self.views = project.state.views;
        self.export_width = project.state.width;
        self.export_height = project.state.height;
        self.full_plot = project.state.full_plot;
        self.thumbnails.clear();
        self.texture = None;
        self.areas.clear();
        self.text_regions.clear();
        self.dirty = true;
        let revisions = self.history.revisions();
        self.data_versions.retain(|r, _| revisions.contains(r));
        while self.data_versions.len() > 8 {
            let oldest = *self.data_versions.keys().next().unwrap();
            self.data_versions.remove(&oldest);
        }
        let retained = self
            .data_versions
            .keys()
            .copied()
            .chain([self.data_revision])
            .collect();
        self.history.retain_revisions(&retained);
        self.status="Data changes applied · Undo restores the previous plotted data; Save Project records the new sources".into();
        self.error = false;
        Ok(())
    }
    fn begin_sources(&mut self, ctx: &egui::Context) {
        if !self.commit_edit() {
            return;
        }
        if let Some(origin) = &self.origin {
            let mut origin = origin.clone();
            origin.retain_imports(&self.fig.plots.iter().map(|p| p.id).collect());
            self.sources = Some(sources::Editor::new(
                origin,
                self.snapshot(),
                Some(self.fig.clone()),
                None,
                ctx,
            ));
        } else {
            self.status = self.origin_error.clone();
            self.error = true;
        }
    }
    fn text(&self, target: TextTarget) -> &str {
        match target {
            TextTarget::FigureTitle => &self.fig.title,
            TextTarget::FigureCaption => &self.fig.caption,
            TextTarget::PanelTitle(i) => &self.fig.plots[i].title,
            TextTarget::PanelCaption(i) => &self.fig.plots[i].caption,
            TextTarget::XLabel(i) => &self.fig.plots[i].xlabel,
            TextTarget::YLabel(i) => &self.fig.plots[i].ylabel,
        }
    }
    fn begin_edit(&mut self, target: TextTarget, position: Option<egui::Pos2>) {
        if self.fig.plots.first().is_some_and(|p| p.qc.margin_off) {
            self.status = "Margin-free strips omit titles, labels, and captions. Reopen without -margin_off to edit text.".into();
            self.error = true;
            return;
        }
        self.editor = Some(TextEditor {
            target,
            value: self.text(target).into(),
            first: true,
            position,
        });
        self.drag_panel = None;
        self.zoom_start = None;
    }
    fn commit_edit(&mut self) -> bool {
        let Some(editor) = &self.editor else {
            return true;
        };
        if editor.value.len() > 16_384 {
            self.status = "Text exceeds 16384 bytes".into();
            self.error = true;
            return false;
        }
        let before = document::Presentation::capture(&self.fig);
        let text = editor.value.clone();
        match editor.target {
            TextTarget::FigureTitle => self.fig.title = text,
            TextTarget::FigureCaption => self.fig.caption = text,
            TextTarget::PanelTitle(i) => self.fig.plots[i].title = text,
            TextTarget::PanelCaption(i) => self.fig.plots[i].caption = text,
            TextTarget::XLabel(i) => self.fig.plots[i].xlabel = text,
            TextTarget::YLabel(i) => self.fig.plots[i].ylabel = text,
        }
        if let Err(e) = self.renderer.scene(
            &self.fig,
            &self.views,
            self.export_width,
            self.export_height,
        ) {
            before.apply(&mut self.fig).expect("existing presentation");
            self.status = format!(
                "Edit would not fit: {e:#}. Shorten the text or cancel and increase image height."
            );
            self.error = true;
            return false;
        }
        self.editor = None;
        self.edited(before);
        self.status = "Text updated · Undo available".into();
        self.error = false;
        true
    }
    fn edit_window(&mut self, ctx: &egui::Context) {
        let Some(editor) = &mut self.editor else {
            return;
        };
        let caption = matches!(
            editor.target,
            TextTarget::FigureCaption | TextTarget::PanelCaption(_)
        );
        let mut done = false;
        let mut cancel = false;
        let mut window = egui::Window::new(if caption {
            "Edit caption"
        } else {
            "Edit title / axis label"
        })
        .id(egui::Id::new("text-editor"))
        .collapsible(false)
        .resizable(true)
        .default_width(400.);
        if let Some(position) = editor.position {
            window = window.default_pos(position);
        }
        window.show(ctx, |ui| {
            let response = if caption {
                ui.add(
                    egui::TextEdit::multiline(&mut editor.value)
                        .desired_width(f32::INFINITY)
                        .desired_rows(5),
                )
            } else {
                ui.add(egui::TextEdit::singleline(&mut editor.value).desired_width(f32::INFINITY))
            };
            if editor.first {
                response.request_focus();
                editor.first = false;
            }
            if !caption && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                done = true;
            }
            if ui.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape)) {
                cancel = true;
            }
            ui.small(if caption {
                "Captions wrap in saved images. Done commits; Escape cancels."
            } else {
                "Enter commits; Escape cancels. Unicode and ^{…} / _{…} supported."
            });
            ui.horizontal(|ui| {
                done |= ui.button("Done").clicked();
                cancel |= ui.button("Cancel").clicked();
            });
        });
        if cancel {
            self.editor = None;
        } else if done {
            self.commit_edit();
        }
    }
    fn project_save(&mut self, save_as: bool) -> bool {
        if !self.commit_style() || !self.commit_edit() {
            return false;
        }
        let Some(origin) = &self.origin else {
            self.status = self.origin_error.clone();
            self.error = true;
            return false;
        };
        let path = if !save_as && self.project_path.is_some() {
            self.project_path.clone()
        } else {
            rfd::FileDialog::new()
                .set_title("Save editable project")
                .add_filter("1dPetePlot project", &["json"])
                .set_file_name(
                    self.project_path
                        .as_ref()
                        .and_then(|p| p.file_name())
                        .and_then(|p| p.to_str())
                        .unwrap_or("Figure.peteplot.json"),
                )
                .save_file()
        };
        let Some(mut path) = path else {
            return false;
        };
        if path.extension().is_none() {
            path.set_extension("peteplot.json");
        }
        let overwrite = path.exists();
        if overwrite
            && rfd::MessageDialog::new()
                .set_title("Replace existing project?")
                .set_description(format!(
                    "Replace {} with this editable figure?",
                    path.display()
                ))
                .set_buttons(rfd::MessageButtons::YesNo)
                .show()
                != rfd::MessageDialogResult::Yes
        {
            return false;
        }
        let state = self.snapshot();
        let layout = if self.fig.plots.is_empty() {
            Ok(())
        } else {
            self.renderer
                .scene(
                    &self.fig,
                    &self.views,
                    self.export_width,
                    self.export_height,
                )
                .map(|_| ())
        };
        let result = layout.and_then(|_| document::save(&path, origin, &state, overwrite));
        match result {
            Ok(()) => {
                self.saved = state;
                self.saved_data_revision = self.data_revision;
                self.project_path = Some(path.clone());
                self.status = format!("Project saved: {}", path.display());
                self.error = false;
                true
            }
            Err(e) => {
                self.status = format!("Project not saved: {e:#}");
                self.error = true;
                false
            }
        }
    }
    fn may_discard(&mut self) -> bool {
        if !self.modified() {
            return true;
        }
        match rfd::MessageDialog::new().set_title("Unsaved figure changes")
            .set_description("Save the editable project before continuing? Yes saves, No discards changes, Cancel keeps this figure open.")
            .set_buttons(rfd::MessageButtons::YesNoCancel).show() {
            rfd::MessageDialogResult::Yes => self.project_save(false),
            rfd::MessageDialogResult::No => true,
            _ => false,
        }
    }
    fn install_project(&mut self, project: document::Loaded, path: PathBuf) {
        self.import.cancel();
        self.sources = None;
        self.data_versions.clear();
        self.data_revision = 0;
        self.saved_data_revision = 0;
        self.dormant.clear();
        self.dormant_views.clear();
        self.next_id = project
            .figure
            .plots
            .iter()
            .map(|p| p.id + 1)
            .max()
            .unwrap_or(0);
        self.fig = project.figure;
        self.origin = Some(project.origin);
        self.origin_error.clear();
        self.views = project.state.views.clone();
        self.export_width = project.state.width;
        self.export_height = project.state.height;
        self.full_plot = project.state.full_plot;
        self.saved = project.state;
        self.project_path = Some(path);
        self.history = Default::default();
        self.editor = None;
        self.styling = None;
        self.duplicate_pending = None;
        self.focus = None;
        self.edit_panel = 0;
        self.drag_panel = None;
        self.zoom_start = None;
        self.context_panel = None;
        self.pending = None;
        self.thumbnails.clear();
        self.texture = None;
        self.areas.clear();
        self.text_regions.clear();
        self.dirty = true;
        self.status = "Project opened · source fingerprints verified".into();
        self.error = false;
    }
    fn project_open(&mut self, ctx: &egui::Context) {
        let Some(path) = rfd::FileDialog::new()
            .set_title("Open editable project")
            .add_filter("1dPetePlot project", &["json"])
            .pick_file()
        else {
            return;
        };
        let result = document::load(&path).and_then(|p| {
            if !p.figure.plots.is_empty() {
                self.renderer
                    .scene(&p.figure, &p.state.views, p.state.width, p.state.height)?;
            }
            Ok(p)
        });
        match result {
            Ok(project) => {
                if self.may_discard() {
                    self.install_project(project, path);
                }
            }
            Err(e) => {
                if let Ok((origin, state)) = document::read_project(&path) {
                    let mut editor = sources::Editor::new(origin, state, None, Some(path), ctx);
                    editor.reject(format!("Project could not be rebuilt: {e:#}"));
                    self.sources = Some(editor);
                }
                self.status = format!("Project not opened: {e:#}");
                self.error = true;
            }
        }
    }

    fn begin_attachments(&mut self, index: usize, ctx: &egui::Context) {
        if !self.commit_edit() {
            return;
        }
        let result = (|| -> Result<()> {
            let plot = self
                .fig
                .plots
                .get(index)
                .ok_or_else(|| anyhow!("select a graph"))?;
            let (o, roles) = self
                .origin
                .as_ref()
                .ok_or_else(|| anyhow!("graph has no reusable source"))?
                .editable(plot.id)?;
            let panel = document::Presentation::capture(&self.fig).panels[index].clone();
            self.import.edit(plot.id, o, roles, panel, ctx)
        })();
        if let Err(e) = result {
            self.status = format!("Cannot edit attachments: {e:#}");
            self.error = true;
        }
    }
    fn replace_import(&mut self, id: usize, prepared: crate::import::Prepared) -> Result<()> {
        prepared.origin.verify_all()?;
        let mut figure = self.fig.clone();
        let mut views = self.views.clone();
        let index = figure
            .plots
            .iter()
            .position(|p| p.id == id)
            .ok_or_else(|| anyhow!("edited graph no longer exists"))?;
        anyhow::ensure!(
            prepared.figure.plots.len() == 1,
            "graph edit needs one graph"
        );
        let mut plot = prepared.figure.plots[0].clone();
        plot.id = id;
        views[index] = plot.bounds;
        figure.plots[index] = plot;
        let mut origin = self
            .origin
            .clone()
            .ok_or_else(|| anyhow!("missing source recipe"))?;
        origin.retain_imports(&self.fig.plots.iter().map(|p| p.id).collect());
        origin.replace_graph(id, prepared.origin, prepared.roles);
        let mut state = self.snapshot();
        state.presentation = document::Presentation::capture(&figure);
        state.views = views;
        state.height = state
            .height
            .max(crate::rsa::minimum_height(&figure, state.width));
        self.apply_data(document::Loaded {
            figure,
            origin,
            state,
        })
    }
    fn open_data(&mut self, ctx: &egui::Context) {
        if self.editor.is_some() || self.styling.is_some() {
            self.status = "Finish editing before importing data".into();
            return;
        }
        if let Some(paths) = rfd::FileDialog::new()
            .set_title("Open data for preview")
            .add_filter("AFNI tables", &["1D", "1d", "tsv"])
            .pick_files()
        {
            self.import.enqueue(paths, ctx);
        }
    }
    fn add_import(&mut self, prepared: crate::import::Prepared) -> Result<()> {
        anyhow::ensure!(
            self.fig.plots.len() + self.dormant.len() < 64,
            "This session retains 64 graphs including Undo history. Save and reopen the project before importing more."
        );
        anyhow::ensure!(
            self.fig.plots.iter().all(|p| !p.qc.margin_off),
            "A margin-free QC strip cannot share a layout. Open a new window to import these data."
        );
        anyhow::ensure!(
            self.origin.is_some(),
            "This figure uses stdin. Save its input to a file and reopen it before importing additional graphs."
        );
        let before = document::Presentation::capture(&self.fig);
        let mut figure = self.fig.clone();
        let mut imported = prepared.figure;
        anyhow::ensure!(
            imported.plots.len() == 1,
            "import preview must contain one graph"
        );
        let mut plot = imported.plots.remove(0);
        plot.id = self.next_id;
        let mut views = self.views.clone();
        views.push(plot.bounds);
        figure.plots.push(plot);
        figure.subtitle = "Imported data · independent graph recipes".into();
        let minimum = ((150. + figure.plots.len().div_ceil(figure.columns) as f64 * 290.)
            * self.export_width as f64
            / 1200.)
            .ceil() as u32;
        let max_height = 8192.min(32_000_000 / self.export_width);
        let mut height = self.export_height.max(minimum);
        loop {
            anyhow::ensure!(
                height <= max_height,
                "The combined layout exceeds image limits. Increase Columns or reduce export width, then retry."
            );
            match self
                .renderer
                .scene(&figure, &views, self.export_width, height)
            {
                Ok(_) => break,
                Err(e) if height + 256 > max_height => return Err(e),
                Err(_) => height += 256,
            }
        }
        self.origin
            .as_mut()
            .unwrap()
            .add(self.next_id, prepared.origin, prepared.roles);
        self.next_id += 1;
        self.fig = figure;
        self.views = views;
        self.export_height = height;
        self.focus = Some(self.fig.plots.len() - 1);
        self.edit_panel = self.fig.plots.len() - 1;
        self.edited(before);
        self.status =
            "Graph added · Undo available · Save Project retains its source and column choices"
                .into();
        self.error = false;
        Ok(())
    }
    fn example(&mut self) {
        if !self.fig.plots.is_empty() {
            return;
        }
        self.fig = crate::model::demo();
        self.views = self.fig.plots.iter().map(|p| p.bounds).collect();
        self.origin = Some(
            document::Origin::capture(&crate::cli::parse(&["-demo".into()]).unwrap()).unwrap(),
        );
        self.next_id = self.fig.plots.len();
        self.data_versions.clear();
        self.data_revision = 0;
        self.saved_data_revision = 0;
        self.history = Default::default();
        self.dormant.clear();
        self.dormant_views.clear();
        self.saved = self.snapshot();
        self.dirty = true;
        self.status = "Synthetic gallery · Open Data adds real data as new graphs".into();
    }

    fn selected_figure(&self, panel: Option<usize>) -> (Figure, Vec<Bounds>) {
        if let Some(i) = panel {
            let mut f = self.fig.clone();
            f.plots = vec![f.plots[i].clone()];
            f.title = f.plots[0].title.clone();
            f.caption.clear();
            f.columns = 1;
            (f, vec![self.views[i]])
        } else {
            (self.fig.clone(), self.views.clone())
        }
    }
    fn displayed(&self) -> (Figure, Vec<Bounds>) {
        self.selected_figure(self.focus)
    }
    fn move_panel(&mut self, from: usize, to: usize) {
        if from == to {
            return;
        }
        let before = document::Presentation::capture(&self.fig);
        let plot = self.fig.plots.remove(from);
        self.fig.plots.insert(to, plot);
        let view = self.views.remove(from);
        self.views.insert(to, view);
        if let Some(i) = self.focus {
            self.focus = Some(if i == from {
                to
            } else if from < i && i <= to {
                i - 1
            } else if to <= i && i < from {
                i + 1
            } else {
                i
            });
        }
        self.edited(before);
        self.edit_panel = to;
        self.context_panel = None;
        self.thumbnails.clear();
        self.dirty = true;
    }
    fn graph_context_menu(&mut self, ui: &mut egui::Ui, panel: Option<usize>) {
        if let Some(i) = panel {
            ui.label(&self.fig.plots[i].title);
            if ui.button("Style this graph…").clicked() {
                self.begin_style(i);
                ui.close();
            }
            if ui.button("Duplicate this graph").clicked() {
                self.duplicate_pending = Some(self.fig.plots[i].id);
                ui.close();
            }
            if ui.button("Edit scientific attachments…").clicked() {
                self.begin_attachments(i, ui.ctx());
                ui.close();
            }
            ui.separator();
            if ui.button("Copy this graph").clicked() {
                self.copy_image(Some(i));
                ui.close();
            }
            if ui.button("Save this graph…").clicked() {
                let focus = self.focus;
                self.focus = Some(i);
                self.request_save();
                self.focus = focus;
                ui.close();
            }
            ui.separator();
        }
        if ui.button("Copy all graphs").clicked() {
            self.copy_image(None);
            ui.close();
        }
        if ui.button("Save all graphs…").clicked() {
            let focus = self.focus;
            self.focus = None;
            self.request_save();
            self.focus = focus;
            ui.close();
        }
        ui.small(if self.full_plot {
            "Full ranges · export dimensions"
        } else {
            "Current view · export dimensions"
        });
    }
    fn sidebar(&mut self, ui: &mut egui::Ui) {
        let titles: Vec<_> = self.fig.plots.iter().map(|p| p.title.clone()).collect();
        if self.thumbnails.len() != titles.len() || titles != self.thumbnail_titles {
            self.thumbnails.clear();
            for i in 0..self.fig.plots.len() {
                let (mut f, _) = self.selected_figure(Some(i));
                // The navigator is a compact plot preview; full text stays on the canvas.
                f.plots[0].caption.clear();
                f.plots[0].title = f.plots[0].title.chars().take(80).collect();
                f.title = f.plots[0].title.clone();
                let result = self
                    .renderer
                    .scene(&f, &[f.plots[0].bounds], 400, 300)
                    .and_then(|scene| self.renderer.raster(&scene));
                match result {
                    Ok(pix) => self.thumbnails.push(ui.ctx().load_texture(
                        format!("thumbnail-{i}"),
                        egui::ColorImage::from_rgba_premultiplied([400, 300], pix.data()),
                        egui::TextureOptions::LINEAR,
                    )),
                    Err(e) => {
                        self.status = format!("Rendering thumbnail: {e:#}");
                        self.error = true;
                        break;
                    }
                }
            }
            self.thumbnail_titles = titles;
        }
        ui.label(
            RichText::new("GRAPHS")
                .size(11.)
                .color(Color32::from_rgb(20, 125, 146)),
        );
        let overview = ui.add_sized(
            [ui.available_width(), 32.],
            egui::Button::new(format!("All graphs · {}", self.fig.plots.len()))
                .selected(self.focus.is_none()),
        );
        if overview.clicked() {
            self.focus = None;
            self.dirty = true;
        }
        overview.context_menu(|ui| self.graph_context_menu(ui, None));
        ui.small("Drag thumbnails to reorder");
        ui.separator();
        let mut reorder = None;
        egui::ScrollArea::vertical()
            .id_salt("graph-list")
            .show(ui, |ui| {
                for i in 0..self.thumbnails.len() {
                    let texture = self.thumbnails[i].clone();
                    let selected = self.focus == Some(i);
                    let width = (ui.available_width() - 16.).max(60.);
                    let card = egui::Frame::new()
                        .inner_margin(6.)
                        .corner_radius(6.)
                        .fill(if selected {
                            Color32::from_rgb(223, 240, 244)
                        } else {
                            Color32::WHITE
                        })
                        .stroke(egui::Stroke::new(
                            if selected { 2. } else { 1. },
                            if selected {
                                Color32::from_rgb(20, 125, 146)
                            } else {
                                Color32::from_rgb(220, 228, 232)
                            },
                        ))
                        .show(ui, |ui| {
                            ui.set_width(width);
                            ui.image((texture.id(), egui::vec2(width, width * 0.75)));
                            ui.add(
                                egui::Label::new(
                                    RichText::new(format!(
                                        "{}  {}",
                                        i + 1,
                                        self.fig.plots[i].title
                                    ))
                                    .size(12.),
                                )
                                .wrap(),
                            );
                        });
                    let response = ui
                        .interact(
                            card.response.rect,
                            egui::Id::new(("graph", i)),
                            egui::Sense::click_and_drag(),
                        )
                        .on_hover_text(
                            "Click to view · drag to rearrange · right-click to copy or save",
                        );
                    if response.clicked() {
                        self.focus = Some(i);
                        self.edit_panel = i;
                        self.dirty = true;
                    }
                    response.dnd_set_drag_payload(i);
                    if response.dnd_hover_payload::<usize>().is_some() {
                        ui.painter().rect_stroke(
                            card.response.rect,
                            6.,
                            egui::Stroke::new(2., Color32::from_rgb(20, 125, 146)),
                            egui::StrokeKind::Inside,
                        );
                    }
                    if let Some(from) = response.dnd_release_payload::<usize>() {
                        reorder = Some((*from, i));
                    }
                    response.context_menu(|ui| self.graph_context_menu(ui, Some(i)));
                    ui.add_space(8.);
                }
            });
        if let Some((from, to)) = reorder {
            self.move_panel(from, to);
        }
    }
    fn copy_image(&mut self, panel: Option<usize>) {
        let result = (|| -> Result<()> {
            let (f, mut views) = self.selected_figure(panel);
            if self.full_plot {
                views = f.plots.iter().map(|p| p.bounds).collect();
            }
            anyhow::ensure!(
                self.export_width as u64 * self.export_height as u64 <= 32_000_000,
                "Export exceeds 32 million pixels; reduce dimensions"
            );
            let scene = self
                .renderer
                .scene(&f, &views, self.export_width, self.export_height)?;
            let pix = self.renderer.raster(&scene)?;
            if self.clipboard.is_none() {
                self.clipboard = Some(arboard::Clipboard::new()?);
            }
            // Clipboard APIs expect straight RGBA; tiny-skia stores premultiplied RGBA.
            let rgba: Vec<u8> = pix
                .pixels()
                .iter()
                .flat_map(|p| {
                    let c = p.demultiply();
                    [c.red(), c.green(), c.blue(), c.alpha()]
                })
                .collect();
            self.clipboard
                .as_mut()
                .unwrap()
                .set_image(arboard::ImageData {
                    width: pix.width() as usize,
                    height: pix.height() as usize,
                    bytes: std::borrow::Cow::Owned(rgba),
                })?;
            Ok(())
        })();
        match result {
            Ok(()) => {
                self.status = format!(
                    "Copied {} · {} × {}",
                    if panel.is_some() {
                        "graph"
                    } else {
                        "all graphs"
                    },
                    self.export_width,
                    self.export_height
                );
                self.error = false;
            }
            Err(e) => {
                self.status = format!("Clipboard unavailable: {e:#}. Use Save Image instead.");
                self.error = true;
            }
        }
    }
    fn action(&mut self, action: &str, ctx: &egui::Context) {
        if self.sources.is_some() && !matches!(action, "quit" | "fullscreen") {
            self.status = "Apply or cancel source changes first".into();
            return;
        }
        if self.styling.is_some() && !matches!(action, "quit" | "fullscreen") {
            self.status = "Apply or cancel styling before changing the figure".into();
            return;
        }
        if self.import.visible()
            && (!matches!(action, "data-open" | "quit" | "fullscreen")
                || self.import.edit_target.is_some() && action == "data-open")
        {
            self.status = "Finish or cancel the import before changing the figure".into();
            return;
        }
        if self.fig.plots.is_empty()
            && matches!(
                action,
                "save" | "copy" | "copy-all" | "reset" | "style" | "duplicate"
            )
        {
            self.status = "Open or drop data to add your first graph".into();
            return;
        }
        if self.editor.is_some() && matches!(action, "undo" | "redo") {
            ctx.input_mut(|i| {
                i.events.push(egui::Event::Key {
                    key: egui::Key::Z,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: if action == "redo" {
                        egui::Modifiers::COMMAND | egui::Modifiers::SHIFT
                    } else {
                        egui::Modifiers::COMMAND
                    },
                })
            });
            return;
        }
        if self.editor.is_some()
            && matches!(
                action,
                "reset" | "overview" | "copy" | "copy-all" | "save" | "style" | "duplicate"
            )
        {
            return;
        }
        match action {
            "style" => self.begin_style(self.focus.unwrap_or(self.edit_panel)),
            "duplicate" => {
                self.duplicate_pending = self
                    .fig
                    .plots
                    .get(self.focus.unwrap_or(self.edit_panel))
                    .map(|p| p.id)
            }
            "data-open" => self.open_data(ctx),
            "project-open" => self.project_open(ctx),
            "sources" => self.begin_sources(ctx),
            "attachments" => self.begin_attachments(self.focus.unwrap_or(self.edit_panel), ctx),
            "project-save" => {
                self.project_save(false);
            }
            "project-save-as" => {
                self.project_save(true);
            }
            "undo" => self.undo(false),
            "redo" => self.undo(true),
            "quit" => {
                if self.may_discard() {
                    self.closing = true;
                    ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                }
            }
            "save" => self.request_save(),
            "copy" => self.copy_image(Some(self.focus.unwrap_or(self.edit_panel))),
            "copy-all" => self.copy_image(None),
            "reset" => self.reset(),
            "fullscreen" => {
                ctx.send_viewport_cmd(egui::ViewportCommand::Fullscreen(!is_fullscreen(ctx)))
            }
            "overview" => {
                self.focus = None;
                self.dirty = true;
            }
            _ => {}
        }
    }
    fn reset(&mut self) {
        if let Some(i) = self.focus {
            self.views[i] = self.fig.plots[i].bounds;
        } else {
            self.views = self.fig.plots.iter().map(|p| p.bounds).collect();
        }
        self.dirty = true;
        self.status = "View reset".into();
        self.error = false;
    }
    fn request_save(&mut self) {
        let initial = PathBuf::from(&self.save_path);
        let ext = if initial
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("svg"))
        {
            "svg"
        } else {
            "png"
        };
        let mut dialog = rfd::FileDialog::new()
            .set_title("Save Plot Image")
            .add_filter(
                if ext == "png" {
                    "PNG image"
                } else {
                    "SVG vector image"
                },
                &[ext],
            )
            .set_file_name(initial.file_name().unwrap_or_default().to_string_lossy());
        if let Some(parent) = initial.parent().filter(|p| !p.as_os_str().is_empty()) {
            dialog = dialog.set_directory(parent);
        }
        if let Some(mut path) = dialog.save_file() {
            if path.extension().is_none() {
                path.set_extension(ext);
            }
            self.save_path = path.to_string_lossy().into_owned();
            self.save_to_path(path);
        } else {
            self.status = "Save cancelled".into();
            self.error = false;
        }
    }
    fn save_to_path(&mut self, path: PathBuf) {
        let (f, mut views) = self.displayed();
        if self.full_plot {
            views = f.plots.iter().map(|p| p.bounds).collect();
        }
        let req = SaveRequest {
            fig: f,
            views,
            width: self.export_width,
            height: self.export_height,
            path,
        };
        if req.path.exists() {
            self.pending = Some(req)
        } else {
            self.save(req, false)
        }
    }
    fn save(&mut self, req: SaveRequest, overwrite: bool) {
        if req.width as u64 * req.height as u64 > 32_000_000 {
            self.status = "Export exceeds 32 million pixels; reduce dimensions".into();
            self.error = true;
            return;
        }
        match self.renderer.export(
            &req.fig, &req.views, req.width, req.height, &req.path, overwrite,
        ) {
            Ok(()) => {
                self.status = format!(
                    "Saved {} × {} · {}",
                    req.width,
                    req.height,
                    req.path.display()
                );
                self.error = false;
            }
            Err(e) => {
                self.status = format!("{e:#}");
                self.error = true;
            }
        }
    }
    fn change_view(&mut self, index: usize, b: Bounds) {
        let orig = self.fig.plots[index].bounds;
        let ratio = (b.x[1] - b.x[0]) / (orig.x[1] - orig.x[0]);
        if b.valid() && (1e-6..1e6).contains(&ratio) {
            self.views[index] = b;
            self.dirty = true;
        }
    }
    fn canvas(&mut self, ui: &mut egui::Ui) {
        if self.fig.plots.is_empty() {
            ui.vertical_centered(|ui| {
                ui.add_space((ui.available_height() * 0.2).max(20.));
                ui.heading(RichText::new("Start with your data").size(30.));
                ui.add_space(12.);
                ui.label("Drop .1D or TSV files here, or choose Open Data.");
                ui.label("Inspect columns, preview a graph, then add it to your figure.");
                ui.add_space(18.);
                if ui
                    .add_sized([230., 40.], egui::Button::new("Open Data…"))
                    .clicked()
                {
                    self.open_data(ui.ctx());
                }
                if ui.button("Explore the synthetic gallery").clicked() {
                    self.example();
                }
            });
            return;
        }
        let available = ui.available_size();
        let ratio = self.export_width as f32 / self.export_height as f32;
        let size = egui::vec2(
            available.x.min(available.y * ratio).max(1.),
            available.y.min(available.x / ratio).max(1.),
        );
        let (outer, _) = ui.allocate_exact_size(available, egui::Sense::hover());
        let rect = egui::Rect::from_center_size(outer.center(), size);
        let response = ui.interact(
            rect,
            egui::Id::new("plot_canvas"),
            egui::Sense::click_and_drag(),
        );
        let max_texture = ui.input(|i| i.max_texture_side).min(4096) as f32;
        let scale = ui
            .ctx()
            .pixels_per_point()
            .min(2.)
            .min(max_texture / size.x.max(size.y));
        let resolution = [
            (size.x * scale).round().max(1.) as u32,
            (size.y * scale).round().max(1.) as u32,
        ];
        if self.dirty || resolution != self.last_size {
            let (fig, views) = self.displayed();
            match self
                .renderer
                .scene(&fig, &views, resolution[0], resolution[1])
                .and_then(|scene| {
                    let pix = self.renderer.raster(&scene)?;
                    Ok((scene, pix))
                }) {
                Ok((scene, pix)) => {
                    let image = egui::ColorImage::from_rgba_premultiplied(
                        [pix.width() as usize, pix.height() as usize],
                        pix.data(),
                    );
                    if let Some(t) = &mut self.texture {
                        t.set(image, egui::TextureOptions::LINEAR)
                    } else {
                        self.texture = Some(ui.ctx().load_texture(
                            "plot",
                            image,
                            egui::TextureOptions::LINEAR,
                        ));
                    }
                    self.text_regions = scene.text_regions;
                    self.areas = scene.areas;
                    self.logical_height = scene.logical_height;
                    self.last_size = resolution;
                    self.dirty = false;
                }
                Err(e) => {
                    self.status = format!("{e:#}");
                    self.error = true;
                    self.texture = None;
                    self.areas.clear();
                    self.text_regions.clear();
                    self.last_size = resolution;
                    self.dirty = false;
                }
            }
        }
        if let Some(t) = &self.texture {
            ui.painter().image(
                t.id(),
                rect,
                egui::Rect::from_min_max(egui::pos2(0., 0.), egui::pos2(1., 1.)),
                Color32::WHITE,
            );
        }
        if self.editor.is_some() || self.styling.is_some() {
            return;
        }
        if self.drag_panel.is_none()
            && let Some(pos) = response.hover_pos()
        {
            let x = (pos.x - rect.left()) as f64 / rect.width() as f64 * 1200.;
            let y = (pos.y - rect.top()) as f64 / rect.height() as f64 * self.logical_height;
            let region = self.text_regions.iter().find(|r| {
                x >= r.rect[0]
                    && x <= r.rect[0] + r.rect[2]
                    && y >= r.rect[1]
                    && y <= r.rect[1] + r.rect[3]
            });
            if let Some(region) = region {
                let target = match (self.focus, region.target) {
                    (Some(i), TextTarget::FigureTitle | TextTarget::PanelTitle(_)) => {
                        TextTarget::PanelTitle(i)
                    }
                    (Some(i), TextTarget::FigureCaption | TextTarget::PanelCaption(_)) => {
                        TextTarget::PanelCaption(i)
                    }
                    (Some(i), TextTarget::XLabel(_)) => TextTarget::XLabel(i),
                    (Some(i), TextTarget::YLabel(_)) => TextTarget::YLabel(i),
                    (_, t) => t,
                };
                ui.ctx().set_cursor_icon(egui::CursorIcon::Text);
                self.hover = if self.text(target).is_empty() {
                    "Click to add a caption or label"
                } else {
                    "Click to edit text"
                }
                .into();
                if response.clicked() {
                    self.begin_edit(target, Some(pos));
                }
                return;
            }
        }
        let mut hit = None;
        self.hover.clear();
        if let Some(pos) = response.hover_pos() {
            let x = (pos.x - rect.left()) as f64 / rect.width() as f64 * 1200.;
            let y = (pos.y - rect.top()) as f64 / rect.height() as f64 * self.logical_height;
            hit = self
                .areas
                .iter()
                .position(|a| a.contains(x, y))
                .map(|i| (i, x, y, pos));
        }
        // Latch the gesture at press, so changing modifiers during a drag is harmless.
        if response.is_pointer_button_down_on() && ui.input(|i| i.pointer.primary_pressed()) {
            self.drag_panel = hit.map(|v| v.0);
            self.zoom_start = if ui.input(|i| i.modifiers.shift) {
                hit.map(|v| v.3)
            } else {
                None
            };
        }
        if response.secondary_clicked() {
            self.context_panel = hit.map(|v| self.focus.unwrap_or(v.0));
        }
        response.context_menu(|ui| self.graph_context_menu(ui, self.context_panel));
        if let Some((i, x, y, pos)) = hit {
            let a = self.areas[i];
            let point = a.data(x, y);
            let index = self.focus.unwrap_or(i);
            self.hover = self.fig.plots[index].hover(point[0], point[1]);
            ui.ctx().set_cursor_icon(if response.dragged() {
                egui::CursorIcon::Grabbing
            } else {
                egui::CursorIcon::Crosshair
            });
            let plot_rect = egui::Rect::from_min_size(
                egui::pos2(
                    rect.left() + (a.x / 1200.) as f32 * rect.width(),
                    rect.top() + (a.y / self.logical_height) as f32 * rect.height(),
                ),
                egui::vec2(
                    (a.w / 1200.) as f32 * rect.width(),
                    (a.h / self.logical_height) as f32 * rect.height(),
                ),
            );
            let stroke = egui::Stroke::new(0.7, Color32::from_rgba_unmultiplied(50, 80, 95, 95));
            ui.painter().line_segment(
                [
                    egui::pos2(pos.x, plot_rect.top()),
                    egui::pos2(pos.x, plot_rect.bottom()),
                ],
                stroke,
            );
            ui.painter().line_segment(
                [
                    egui::pos2(plot_rect.left(), pos.y),
                    egui::pos2(plot_rect.right(), pos.y),
                ],
                stroke,
            );
            if response.double_clicked() {
                self.views[index] = self.fig.plots[index].bounds;
                self.dirty = true;
            }
        }
        if let (Some(i), Some(start), Some(end)) = (
            self.drag_panel,
            self.zoom_start,
            ui.input(|i| i.pointer.latest_pos()),
        ) {
            let a = self.areas[i];
            let plot_rect = egui::Rect::from_min_max(
                egui::pos2(
                    rect.left() + (a.x / 1200.) as f32 * rect.width(),
                    rect.top() + (a.y / self.logical_height) as f32 * rect.height(),
                ),
                egui::pos2(
                    rect.left() + ((a.x + a.w) / 1200.) as f32 * rect.width(),
                    rect.top() + ((a.y + a.h) / self.logical_height) as f32 * rect.height(),
                ),
            );
            let selection = egui::Rect::from_two_pos(plot_rect.clamp(start), plot_rect.clamp(end));
            ui.painter().rect(
                selection,
                0.,
                Color32::from_rgba_unmultiplied(20, 125, 146, 35),
                egui::Stroke::new(1.5, Color32::from_rgb(20, 125, 146)),
                egui::StrokeKind::Inside,
            );
            if response.drag_stopped_by(egui::PointerButton::Primary)
                && selection.width() > 6.
                && selection.height() > 6.
            {
                let data = |p: egui::Pos2| {
                    a.data(
                        (p.x - rect.left()) as f64 / rect.width() as f64 * 1200.,
                        (p.y - rect.top()) as f64 / rect.height() as f64 * self.logical_height,
                    )
                };
                let lo = data(selection.left_bottom());
                let hi = data(selection.right_top());
                self.change_view(
                    self.focus.unwrap_or(i),
                    Bounds {
                        x: [lo[0], hi[0]],
                        y: [lo[1], hi[1]],
                    },
                );
            }
        }
        if response.dragged_by(egui::PointerButton::Primary)
            && self.zoom_start.is_none()
            && let Some(i) = self.drag_panel
        {
            let a = self.areas[i];
            let index = self.focus.unwrap_or(i);
            let delta = response.drag_delta();
            let dx =
                -(delta.x / rect.width()) as f64 * 1200. / a.w * (a.bounds.x[1] - a.bounds.x[0]);
            let dy = (delta.y / rect.height()) as f64 * self.logical_height / a.h
                * (a.bounds.y[1] - a.bounds.y[0]);
            self.change_view(index, self.views[index].pan(dx, dy));
        }
        if ui.input(|i| i.pointer.primary_released()) {
            self.drag_panel = None;
            self.zoom_start = None;
        }
    }
}
impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = ui.ctx().clone();
        let was_importing = self.import.visible();
        if was_importing
            && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape))
        {
            self.import.cancel();
        }
        if self.pending.is_none()
            && self.editor.is_none()
            && self.styling.is_none()
            && !was_importing
            && self.sources.is_none()
        {
            exit_fullscreen_on_escape(&ctx);
        }
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if !dropped.is_empty() {
            if self.editor.is_some()
                || self.styling.is_some()
                || self.pending.is_some()
                || self.sources.is_some()
                || self.import.edit_target.is_some()
            {
                self.status =
                    "Finish the current edit/save dialog, then drop the file again".into();
            } else {
                let paths: Vec<_> = dropped
                    .into_iter()
                    .map(|f| f.path().to_path_buf())
                    .collect();
                if paths.is_empty() {
                    self.status = "Drop local .1D or TSV files with filesystem paths".into();
                    self.error = true;
                } else {
                    self.import.enqueue(paths, &ctx);
                }
            }
        }
        if ctx.input(|i| i.viewport().close_requested()) && !self.closing {
            ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
            if self.may_discard() {
                self.closing = true;
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
        }
        let name = self
            .project_path
            .as_ref()
            .and_then(|p| p.file_name())
            .and_then(|p| p.to_str())
            .unwrap_or("Untitled");
        let title = format!(
            "{}{} — 1dPetePlot",
            name,
            if self.modified() { " •" } else { "" }
        );
        if title != self.window_title {
            ctx.send_viewport_cmd(egui::ViewportCommand::Title(title.clone()));
            self.window_title = title;
        }
        if self.editor.is_some() {
            ui.disable();
        }
        #[cfg(target_os = "macos")]
        if let Some(menu) = &self.menu {
            menu.set_fullscreen(is_fullscreen(&ctx));
        }
        #[cfg(not(target_os = "macos"))]
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::F11)) {
            self.action("fullscreen", &ctx);
        }
        #[cfg(target_os = "macos")]
        while let Some(action) = self.menu.as_ref().and_then(|m| m.next()) {
            self.action(&action, &ctx);
        }
        #[cfg(not(target_os = "macos"))]
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    for (label, action) in [
                        ("Open Data…    Ctrl+I", "data-open"),
                        ("Open Project…    Ctrl+O", "project-open"),
                        ("Save Project…    Ctrl+Shift+S", "project-save"),
                        ("Save Project As…", "project-save-as"),
                        ("Reload / relink sources…", "sources"),
                    ] {
                        if ui.button(label).clicked() {
                            ui.close();
                            self.action(action, &ctx);
                        }
                    }
                    ui.separator();
                    if ui.button("Save Image…    Ctrl+S").clicked() {
                        ui.close();
                        self.action("save", &ctx);
                    }
                });
                ui.menu_button("Edit", |ui| {
                    for (label, action) in [
                        ("Undo    Ctrl+Z", "undo"),
                        ("Redo    Ctrl+Shift+Z", "redo"),
                        ("Style graph…", "style"),
                        ("Edit scientific attachments…", "attachments"),
                        ("Duplicate graph    Ctrl+D", "duplicate"),
                    ] {
                        if ui.button(label).clicked() {
                            ui.close();
                            self.action(action, &ctx);
                        }
                    }
                    ui.separator();
                    if ui.button("Copy selected graph").clicked() {
                        ui.close();
                        self.action("copy", &ctx);
                    }
                    if ui.button("Copy all graphs").clicked() {
                        ui.close();
                        self.action("copy-all", &ctx);
                    }
                });
                ui.menu_button("View", |ui| {
                    if ui
                        .button(if is_fullscreen(&ctx) {
                            "Exit full screen   F11"
                        } else {
                            "Enter full screen   F11"
                        })
                        .clicked()
                    {
                        ui.close();
                        self.action("fullscreen", &ctx);
                    }
                    if ui.button("Reset view").clicked() {
                        ui.close();
                        self.action("reset", &ctx);
                    }
                    if ui.button("Show all graphs").clicked() {
                        ui.close();
                        self.action("overview", &ctx);
                    }
                });
            });
        });
        if self.styling.is_some() {
            self.style_workspace(ui);
            return;
        }
        if let Some(editor) = &mut self.sources {
            let recovery = editor.recovery.clone();
            let (cancel, project) = editor.show(ui);
            if cancel {
                self.sources = None;
            }
            if let Some(project) = project {
                if let Some(path) = recovery {
                    if self.may_discard() {
                        self.install_project(project, path);
                        self.saved_data_revision = u64::MAX;
                        self.status="Recovered project preview applied · Save Project to retain the new source references".into();
                    }
                } else {
                    match self.apply_data(project) {
                        Ok(()) => self.sources = None,
                        Err(e) => {
                            self.status = format!("Apply failed: {e:#}");
                            if let Some(editor) = &mut self.sources {
                                editor.reject(self.status.clone());
                            }
                            self.error = true;
                        }
                    }
                }
            }
            return;
        }
        // Import replaces the workspace so its approval controls use the whole window.
        if self.import.visible() {
            if let Some(prepared) = self.import.show(ui) {
                let edit = self.import.edit_target;
                let result = if let Some(id) = edit {
                    self.replace_import(id, *prepared)
                } else {
                    self.add_import(*prepared)
                };
                match result {
                    Ok(()) => if edit.is_some(){self.import.cancel()}else{self.import.advance()},
                    Err(e) => self.import.reject(format!(
                        "Graph was not added: {e:#}. Cancel import to adjust the layout, then retry."
                    )),
                }
                ctx.request_repaint();
            }
            return;
        }
        egui::Panel::bottom("status")
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(236, 242, 244))
                    .inner_margin(10.),
            )
            .show(ui, |ui| {
                ui.label(
                    RichText::new(if self.hover.is_empty() {
                        &self.status
                    } else {
                        &self.hover
                    })
                    .size(12.)
                    .color(if self.error && self.hover.is_empty() {
                        Color32::from_rgb(175, 50, 50)
                    } else {
                        Color32::from_rgb(57, 79, 94)
                    }),
                );
            });
        egui::Panel::left("graphs")
            .default_size(184.)
            .size_range(150.0..=280.0)
            .resizable(true)
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(237, 242, 244))
                    .inner_margin(12.),
            )
            .show(ui, |ui| self.sidebar(ui));
        egui::Panel::right("inspector").exact_size(252.).frame(egui::Frame::new().fill(Color32::WHITE).inner_margin(18.)).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                if ui.add_sized([210.,32.], egui::Button::new("Open Data…")).clicked() { self.open_data(&ctx); }
                if let Some(error) = &self.import.error && !self.import.visible() { ui.colored_label(Color32::from_rgb(170,55,45), error); }
                ui.separator();
                ui.label(RichText::new(if self.modified() { "PROJECT · Modified" } else { "PROJECT" }).size(11.).color(Color32::from_rgb(20,125,146)));
                ui.horizontal(|ui| {
                    if ui.button("Open…").clicked() { self.project_open(ui.ctx()); }
                    if ui.button("Save Project…").clicked() { self.project_save(false); }
                });
                ui.horizontal(|ui| {
                    if ui.add_enabled(self.history.can_undo(), egui::Button::new("Undo")).clicked() { self.undo(false); }
                    if ui.add_enabled(self.history.can_redo(), egui::Button::new("Redo")).clicked() { self.undo(true); }
                });
                if ui.add_enabled(self.project_path.is_some() && !self.modified(), egui::Button::new("Copy replay command")).on_hover_text("Save the project first; the command exports all panels using saved settings.").clicked()
                    && let Some(path) = self.project_path.as_ref().and_then(|p| std::fs::canonicalize(p).ok()) {
                    let quoted = path.to_string_lossy().replace('\'', "'\"'\"'");
                    ctx.copy_text(format!("1dPetePlot -recipe '{quoted}' -prefix figure.png"));
                    self.status = "Copied command for the saved project".into(); self.error = false;
                }
                if self.fig.plots.is_empty() { ui.small("Your figure is empty. Add graphs by opening or dropping data files."); return; }
                ui.small("Click titles, axis labels, or captions to edit.");
                ui.separator();
                let selected = self.focus.unwrap_or(self.edit_panel);
                ui.label(format!("GRAPH {} · STYLE & DUPLICATE", selected + 1));
                ui.small(&self.fig.plots[selected].title);
                if ui.add_sized([210.,32.], egui::Button::new("Style graph…")).clicked() { self.begin_style(selected); }
                if ui.add_sized([210.,32.], egui::Button::new("Duplicate graph")).clicked() { self.duplicate_pending = Some(self.fig.plots[selected].id); }
                ui.separator();
                ui.label(RichText::new("VIEW").size(11.).color(Color32::from_rgb(20,125,146)));
                if is_fullscreen(&ctx) && ui.add_sized([210., 30.], egui::Button::new("Exit full screen   Esc")).clicked() {
                    ctx.send_viewport_cmd(egui::ViewportCommand::Fullscreen(false));
                }
                if ui.add_sized([210.,30.], egui::Button::new("Reset view   R")).clicked() { self.reset(); }
                ui.label(RichText::new("Drag to pan · Shift-drag to zoom\nDouble-click a graph to reset\nRight-click to copy or save").size(12.).color(Color32::from_rgb(100,117,132)));
                ui.separator();
                ui.collapsing("Titles & arrangement", |ui| {
                    if ui.button("Edit figure title…").clicked() { self.begin_edit(TextTarget::FigureTitle, None); }
                    if ui.button("Edit figure caption…").clicked() { self.begin_edit(TextTarget::FigureCaption, None); }
                    ui.horizontal(|ui| {
                        ui.label("Columns");
                        let before = document::Presentation::capture(&self.fig);
                        if ui.add(egui::DragValue::new(&mut self.fig.columns).range(1..=self.fig.plots.len().min(4))).changed() { self.edited(before); }
                    });
                    egui::ComboBox::from_id_salt("edit-panel").selected_text(format!("Graph {}", self.edit_panel + 1)).show_ui(ui, |ui| {
                        for i in 0..self.fig.plots.len() { ui.selectable_value(&mut self.edit_panel, i, format!("{} · {}", i + 1, self.fig.plots[i].title)); }
                    });
                    for (label, target) in [("Edit graph title…", TextTarget::PanelTitle(self.edit_panel)), ("Edit graph caption…", TextTarget::PanelCaption(self.edit_panel)), ("Edit x label…", TextTarget::XLabel(self.edit_panel)), ("Edit y label…", TextTarget::YLabel(self.edit_panel))] {
                        if ui.button(label).clicked() { self.begin_edit(target, None); }
                    }
                    ui.horizontal(|ui| {
                        let i = self.edit_panel;
                        if ui.add_enabled(i > 0, egui::Button::new("← Earlier")).clicked() { self.move_panel(i, i - 1); }
                        if ui.add_enabled(i + 1 < self.fig.plots.len(), egui::Button::new("Later →")).clicked() { self.move_panel(i, i + 1); }
                    });
                    ui.small("Changes apply to images you save or copy.");
                });
                ui.separator();
                ui.label(RichText::new("SAVE & COPY").size(11.).color(Color32::from_rgb(20,125,146)));
                ui.radio_value(&mut self.full_plot, false, "Current view");
                ui.radio_value(&mut self.full_plot, true, "Full plot");
                ui.label("Dimensions in pixels");
                ui.horizontal(|ui| {
                    let w = ui.add(egui::DragValue::new(&mut self.export_width).range(400..=8192).speed(10.));
                    ui.label("×");
                    let h = ui.add(egui::DragValue::new(&mut self.export_height).range(if self.fig.plots.first().is_some_and(|p| p.qc.margin_off) {32}else{300}..=8192).speed(10.));
                    if w.changed() || h.changed() { self.dirty = true; }
                });
                ui.horizontal(|ui| { for ext in ["png", "svg"] {
                    if ui.selectable_label(self.save_path.to_lowercase().ends_with(&format!(".{ext}")), ext.to_uppercase()).clicked() {
                        let mut p = PathBuf::from(&self.save_path); p.set_extension(ext); self.save_path = p.display().to_string();
                    }
                }});
                if ui.add_sized([210.,34.], egui::Button::new(RichText::new("Save Image…").color(Color32::WHITE)).fill(Color32::from_rgb(20,125,146))).clicked() { self.request_save(); }
                ui.small(if self.focus.is_some() { "Saves the selected graph." } else { "Saves all graphs in this arrangement." });
                ui.label(RichText::new("PNG for sharing. SVG for scalable figures. Clipboard copies use a bitmap at these dimensions.").size(11.).color(Color32::from_rgb(110,126,137)));
                ui.separator();
                let zoomed = self.views.iter().zip(&self.fig.plots).any(|(v,p)| *v != p.bounds);
                ui.label(RichText::new(if zoomed { "VIEW MODIFIED" } else { "FULL DATA RANGE" }).size(10.).color(Color32::from_rgb(20,125,146)));
                if self.fig.synthetic { ui.small("Synthetic examples; not experimental findings."); }
            });
        });
        egui::CentralPanel::default()
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(243, 246, 247))
                    .inner_margin(8.),
            )
            .show(ui, |ui| self.canvas(ui));
        self.edit_window(&ctx);
        if let Some(id) = self.duplicate_pending.take()
            && let Some(index) = self.fig.plots.iter().position(|p| p.id == id)
            && let Err(e) = self.duplicate_graph(index)
        {
            self.status = format!("Graph was not duplicated: {e:#}");
            self.error = true;
        }
        if !ctx.input(|i| i.raw.hovered_files.is_empty()) {
            egui::Area::new(egui::Id::new("drop-hint"))
                .order(egui::Order::Tooltip)
                .anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO)
                .show(&ctx, |ui| {
                    egui::Frame::popup(ui.style()).show(ui, |ui| {
                        ui.heading("Drop to inspect and preview");
                        ui.label("Files are added only after you choose Add graph.");
                    });
                });
        }
        if !self.import.visible()
            && self.editor.is_none()
            && self.styling.is_none()
            && !ctx.egui_wants_keyboard_input()
            && ctx.input(|i| i.key_pressed(egui::Key::R))
        {
            self.reset();
        }
        #[cfg(not(target_os = "macos"))]
        if self.editor.is_none() && self.styling.is_none() {
            for (key, modifiers, action) in [
                (
                    egui::Key::S,
                    egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                    "project-save",
                ),
                (egui::Key::S, egui::Modifiers::COMMAND, "save"),
                (egui::Key::O, egui::Modifiers::COMMAND, "project-open"),
                (egui::Key::I, egui::Modifiers::COMMAND, "data-open"),
                (
                    egui::Key::Z,
                    egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                    "redo",
                ),
                (egui::Key::Z, egui::Modifiers::COMMAND, "undo"),
                (egui::Key::D, egui::Modifiers::COMMAND, "duplicate"),
            ] {
                if ctx.input_mut(|i| i.consume_key(modifiers, key)) {
                    self.action(action, &ctx);
                }
            }
        }
        if let Some(req) = self.pending.clone() {
            egui::Window::new("Replace existing image?")
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO)
                .show(&ctx, |ui| {
                    ui.label(format!("{} already exists.", req.path.display()));
                    ui.label("Replace it with this figure?");
                    ui.horizontal(|ui| {
                        if ui.button("Cancel").clicked() {
                            self.pending = None;
                            self.status = "Save cancelled; existing image preserved".into();
                            self.error = false;
                        }
                        if ui.button("Replace").clicked() {
                            self.pending = None;
                            self.save(req.clone(), true);
                        }
                    });
                });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::App as _;
    fn app() -> App {
        App::new(crate::model::demo(), Options::default(), Renderer::new())
    }
    #[test]
    fn import_preview_add_undo_reorder_and_replay_preserve_pixels_and_sources() {
        let temp = tempfile::tempdir().unwrap();
        let data = temp.path().join("data.tsv");
        std::fs::write(&data, "roi\tx\ty\nA\t0\t1\nB\t1\t3\nC\t2\t2\n").unwrap();
        let table = crate::import::inspect(&data, || false).unwrap();
        let mut a = App::new(Figure::empty(), Options::default(), Renderer::new());
        let mut mapping = crate::import::Mapping::for_kind(&table, crate::import::Kind::Scatter);
        // Numeric source columns 1 and 2; categorical ROI stays unmodified.
        let preview = crate::import::prepare(&table, &mapping, || false).unwrap();
        let preview_data = preview.rgba.clone();
        let preview_bounds = preview.figure.plots[0].bounds;
        a.add_import(preview).unwrap();
        let (fig, views) = a.selected_figure(Some(0));
        assert_eq!(views, [preview_bounds]);
        // Preview figure title/subtitle are presentation; the panel and its samples are identical.
        let mut expected = fig.clone();
        expected.title = a.fig.plots[0].title.clone();
        expected.subtitle = "scatter · 1 input file".into();
        let pixels = a
            .renderer
            .raster(&a.renderer.scene(&expected, &views, 900, 600).unwrap())
            .unwrap();
        assert_eq!(pixels.data(), preview_data);
        let first_id = a.fig.plots[0].id;
        mapping = crate::import::Mapping::for_kind(&table, crate::import::Kind::Histogram);
        a.add_import(crate::import::prepare(&table, &mapping, || false).unwrap())
            .unwrap();
        let second_id = a.fig.plots[1].id;
        assert_ne!(first_id, second_id);
        a.move_panel(1, 0);
        a.undo(false);
        a.undo(false);
        assert_eq!(a.fig.plots.len(), 1);
        assert_eq!(a.fig.plots[0].id, first_id);
        let path = temp.path().join("one.peteplot.json");
        document::save(&path, a.origin.as_ref().unwrap(), &a.snapshot(), false).unwrap();
        assert_eq!(document::load(&path).unwrap().figure.plots.len(), 1);
        a.undo(true);
        assert_eq!(a.fig.plots.len(), 2);
        assert_eq!(a.fig.plots[1].id, second_id);
        a.move_panel(1, 0);
        a.begin_edit(TextTarget::PanelCaption(0), None);
        a.editor.as_mut().unwrap().value = "Imported histogram".into();
        assert!(a.commit_edit());
        let path = temp.path().join("mixed.peteplot.json");
        document::save(&path, a.origin.as_ref().unwrap(), &a.snapshot(), false).unwrap();
        let loaded = document::load(&path).unwrap();
        assert_eq!(loaded.state, a.snapshot());
        let raster = |fig: &Figure| {
            a.renderer
                .raster(
                    &a.renderer
                        .scene(fig, &a.views, a.export_width, a.export_height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(raster(&loaded.figure).data(), raster(&a.fig).data());
        std::fs::write(data, "roi\tx\ty\nA\t0\t8\nB\t1\t3\nC\t2\t2\n").unwrap();
        assert!(document::load(&path).is_err());
    }
    #[test]
    fn rdm_add_undo_and_project_replay_preserve_label_sources_and_pixels() {
        let d = tempfile::tempdir().unwrap();
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/rdm");
        let labels = d.path().join("conditions.txt");
        std::fs::copy(base.join("conditions.txt"), &labels).unwrap();
        let table = crate::import::inspect(&base.join("model.1D"), || false).unwrap();
        let mut mapping = crate::import::Mapping::for_kind(&table, crate::import::Kind::Rdm);
        mapping.rdm.labels = labels.to_str().unwrap().into();
        mapping.rdm.clim = "0:1".into();
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        app.add_import(crate::import::prepare(&table, &mapping, || false).unwrap())
            .unwrap();
        app.undo(false);
        assert!(app.fig.plots.is_empty());
        app.undo(true);
        assert_eq!(app.fig.plots.len(), 1);
        let path = d.path().join("rdm.json");
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        let loaded = document::load(&path).unwrap();
        let pixels = |f: &Figure| {
            app.renderer
                .raster(
                    &app.renderer
                        .scene(f, &app.views, app.export_width, app.export_height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(pixels(&app.fig).data(), pixels(&loaded.figure).data());
        let prepared = crate::import::prepare(&table, &mapping, || false).unwrap();
        std::fs::write(&labels, "Changed\nlabels\n").unwrap();
        assert!(prepared.origin.verify().is_err()); // Add confirmation also checks attached labels.
        assert!(document::load(&path).is_err());
    }
    #[test]
    fn reloaded_data_and_text_edits_undo_redo_in_order_and_reject_stale_apply() {
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("values.1D");
        std::fs::write(&p, "1\n2\n3\n").unwrap();
        let o = crate::cli::parse(&["-infiles".into(), p.to_string_lossy().into(), "-one".into()])
            .unwrap();
        let f = crate::model::from_options(&o).unwrap();
        let mut app = App::new(f, o, Renderer::new());
        let before = document::Presentation::capture(&app.fig);
        app.fig.plots[0].title = "My scientific title".into();
        app.edited(before);
        std::fs::write(&p, "4\n5\n6\n7\n").unwrap();
        let choices = vec![sources::Choice {
            original: p.to_string_lossy().into(),
            replacement: p.to_string_lossy().into(),
            status: "Changed".into(),
        }];
        let preview = sources::prepare(
            app.origin.as_ref().unwrap(),
            &app.snapshot(),
            Some(&app.fig),
            &choices,
            false,
        )
        .unwrap();
        app.apply_data(preview.project).unwrap();
        assert!(app.modified());
        let before = document::Presentation::capture(&app.fig);
        app.fig.plots[0].caption = "After reloading".into();
        app.edited(before);
        app.undo(false);
        assert!(app.fig.plots[0].caption.is_empty());
        app.undo(false);
        let crate::model::Content::Lines(ss) = &app.fig.plots[0].content else {
            panic!()
        };
        assert_eq!(ss[0].points[0][1], 1.);
        assert_eq!(ss[0].points.len(), 3);
        assert_eq!(app.fig.plots[0].title, "My scientific title");
        assert!(app.origin.as_ref().unwrap().verify().is_err());
        app.undo(true);
        let crate::model::Content::Lines(ss) = &app.fig.plots[0].content else {
            panic!()
        };
        assert_eq!(ss[0].points.len(), 4);
        assert_eq!(ss[0].points[0][1], 4.);
        app.undo(true);
        assert_eq!(app.fig.plots[0].caption, "After reloading");
        let preview = sources::prepare(
            app.origin.as_ref().unwrap(),
            &app.snapshot(),
            Some(&app.fig),
            &choices,
            false,
        )
        .unwrap();
        let old = app.snapshot();
        std::fs::write(&p, "9\n").unwrap();
        assert!(app.apply_data(preview.project).is_err());
        assert_eq!(app.snapshot(), old);
    }
    #[test]
    fn attachment_edits_fork_duplicates_and_replay_with_data_undo() {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/attachments");
        let table = crate::import::inspect(&base.join("signals.tsv"), || false).unwrap();
        let mut m = crate::import::Mapping::for_kind(&table, crate::import::Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.uncertainty = "bounds".into();
        m.attachments.bands.insert(1, [Some(2), Some(3)]);
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        app.add_import(crate::import::prepare(&table, &m, || false).unwrap())
            .unwrap();
        app.duplicate_graph(0).unwrap();
        let (o, roles) = app.origin.as_ref().unwrap().editable(0).unwrap();
        let panel = document::Presentation::capture(&app.fig).panels[0].clone();
        let mut edit = crate::import::Mapping::restore(&table, &o, &roles, panel).unwrap();
        assert_eq!(edit.attachments.bands[&1], [Some(2), Some(3)]);
        edit.attachments.uncertainty = "none".into();
        edit.attachments.censor_column = Some(6);
        edit.attachments.xmat = base.join("X.mat.1D").to_string_lossy().into();
        assert!(crate::import::prepare(&table, &edit, || false).is_err());
        edit.attachments.aligned = true;
        app.replace_import(0, crate::import::prepare(&table, &edit, || false).unwrap())
            .unwrap();
        let has_band = |app: &App, id| {
            let crate::model::Content::Lines(ss) =
                &app.fig.plots.iter().find(|p| p.id == id).unwrap().content
            else {
                panic!()
            };
            ss[0].band.is_some()
        };
        assert!(!has_band(&app, 0));
        assert!(has_band(&app, 1));
        assert!(!app.fig.plots[0].qc.annotations.is_empty());
        app.duplicate_graph(0).unwrap();
        assert!(!has_band(&app, 2));
        app.undo(false);
        app.undo(false);
        assert!(has_band(&app, 0));
        app.undo(true);
        assert!(!has_band(&app, 0));
        app.undo(true);
        assert!(!has_band(&app, 2));
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("edited.json");
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        let loaded = document::load(&path).unwrap();
        let r = Renderer::new();
        let pix = |f: &Figure| {
            r.raster(
                &r.scene(f, &app.views, app.export_width, app.export_height)
                    .unwrap(),
            )
            .unwrap()
        };
        assert_eq!(pix(&app.fig).data(), pix(&loaded.figure).data());
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&std::fs::read_to_string(&path).unwrap())
                .unwrap()["version"],
            4
        );
        // Remove the design from both edited copies. It must stop being a replay
        // dependency, while the unchanged original copy still retains its bands.
        for id in [0, 2] {
            let (o, roles) = app.origin.as_ref().unwrap().editable(id).unwrap();
            let panel = document::Presentation::capture(&app.fig)
                .panels
                .into_iter()
                .find(|p| p.id == id)
                .unwrap();
            let mut edit = crate::import::Mapping::restore(&table, &o, &roles, panel).unwrap();
            edit.attachments.xmat.clear();
            app.replace_import(id, crate::import::prepare(&table, &edit, || false).unwrap())
                .unwrap();
        }
        assert!(
            !app.origin
                .as_ref()
                .unwrap()
                .source_paths()
                .iter()
                .any(|p| p.ends_with("X.mat.1D"))
        );
        assert!(has_band(&app, 1));
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), true).unwrap();
        assert!(document::load(&path).is_ok());
    }
    #[test]
    fn long_rsa_imports_add_undo_and_replay_keep_categories_and_bands() {
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        for (name, kind) in [
            ("series", crate::import::Kind::Line),
            ("runs-joint", crate::import::Kind::Scatter),
        ] {
            let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join(format!("examples/rsa/{name}.rsa.1D"));
            let table = crate::import::inspect(&p, || false).unwrap();
            let mut mapping = crate::import::Mapping::for_kind(&table, kind);
            mapping.roi_filter = "12".into();
            mapping.rsa_band = true;
            app.add_import(crate::import::prepare(&table, &mapping, || false).unwrap())
                .unwrap();
        }
        app.undo(false);
        assert_eq!(app.fig.plots.len(), 1);
        app.undo(true);
        assert_eq!(app.fig.plots.len(), 2);
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("long.json");
        document::save(&p, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        let loaded = document::load(&p).unwrap();
        let pixels = |f: &Figure| {
            app.renderer
                .raster(
                    &app.renderer
                        .scene(f, &app.views, app.export_width, app.export_height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(pixels(&app.fig).data(), pixels(&loaded.figure).data());
    }
    #[test]
    fn rsa_imports_add_undo_and_replay_their_semantic_recipes() {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/rsa/ordinary.rsa.1D");
        let table = crate::import::inspect(&path, || false).unwrap();
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        for kind in [
            crate::import::Kind::Histogram,
            crate::import::Kind::Scatter,
            crate::import::Kind::Heat,
            crate::import::Kind::Intervals,
        ] {
            let mut mapping = crate::import::Mapping::for_kind(&table, kind);
            mapping.roi_filter = "3,12,27".into();
            let preview = crate::import::prepare(&table, &mapping, || false).unwrap();
            app.add_import(preview).unwrap();
        }
        assert_eq!(app.fig.plots.len(), 4);
        app.undo(false);
        assert_eq!(app.fig.plots.len(), 3);
        app.undo(true);
        assert_eq!(app.fig.plots.len(), 4);
        let d = tempfile::tempdir().unwrap();
        let project = d.path().join("rsa.json");
        document::save(
            &project,
            app.origin.as_ref().unwrap(),
            &app.snapshot(),
            false,
        )
        .unwrap();
        let loaded = document::load(&project).unwrap();
        assert_eq!(loaded.figure.plots.len(), 4);
        assert!(loaded.figure.plots[3].hover(0.68, 2.5).contains("95% CI"));
        assert!(loaded.figure.plots[1].hover(0.68, 0.21).contains("ROI 3"));
        assert_eq!(
            document::Presentation::capture(&loaded.figure),
            document::Presentation::capture(&app.fig)
        );
    }
    #[test]
    fn empty_workspace_and_dropped_files_do_not_mutate_until_add() {
        #[derive(Debug)]
        struct DropPath(PathBuf);
        impl egui::DroppedFile for DropPath {
            fn path(&self) -> &std::path::Path {
                &self.0
            }
            fn bytes(&self) -> std::result::Result<Vec<u8>, String> {
                std::fs::read(&self.0).map_err(|e| e.to_string())
            }
        }
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("rows.1D");
        std::fs::write(&path, "0\n1\n2\n").unwrap();
        let ctx = egui::Context::default();
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        let before = app.snapshot();
        let mut frame = eframe::Frame::_new_kittest();
        let mut input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(850., 600.),
            )),
            ..Default::default()
        };
        input
            .dropped_files
            .push(std::sync::Arc::new(DropPath(path)));
        ctx.run_ui(input, |ui| app.ui(ui, &mut frame))
            .drop_without_applying_deltas();
        assert!(app.import.visible());
        assert_eq!(app.snapshot(), before);
        let page = app.import.page_rect.unwrap();
        assert_eq!(page.left(), 0.);
        assert_eq!(page.right(), 850.);
        assert_eq!(page.bottom(), 600.);
        assert!(page.top() < 40.); // Only the Linux menu may occupy space above import.
        ctx.run_ui(
            egui::RawInput {
                events: vec![egui::Event::Key {
                    key: egui::Key::Escape,
                    physical_key: None,
                    pressed: true,
                    repeat: false,
                    modifiers: egui::Modifiers::NONE,
                }],
                ..Default::default()
            },
            |ui| app.ui(ui, &mut frame),
        )
        .drop_without_applying_deltas();
        assert!(!app.import.visible());
        assert_eq!(app.snapshot(), before);
        let save = temp.path().join("empty.json");
        document::save(&save, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        assert!(document::load(&save).unwrap().figure.plots.is_empty());
        for action in ["save", "copy", "copy-all", "reset"] {
            app.action(action, &ctx);
        }
        assert!(app.fig.plots.is_empty());
    }
    #[test]
    fn undo_first_import_restores_empty_document_and_redo_branch_releases_data() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("sample.1D");
        std::fs::write(&path, "1\n2\n3\n").unwrap();
        let table = crate::import::inspect(&path, || false).unwrap();
        let mut a = App::new(Figure::empty(), Options::default(), Renderer::new());
        let start = document::Presentation::capture(&a.fig);
        let mapping = crate::import::Mapping::for_kind(&table, crate::import::Kind::Line);
        a.add_import(crate::import::prepare(&table, &mapping, || false).unwrap())
            .unwrap();
        let first_id = a.fig.plots[0].id;
        a.undo(false);
        assert_eq!(document::Presentation::capture(&a.fig), start);
        assert!(a.dormant.contains_key(&first_id));
        a.undo(true);
        assert_eq!(a.fig.plots[0].id, first_id);
        a.undo(false);
        a.add_import(crate::import::prepare(&table, &mapping, || false).unwrap())
            .unwrap();
        assert_ne!(a.fig.plots[0].id, first_id);
        assert!(!a.history.can_redo());
        assert!(a.dormant.is_empty());
    }
    #[test]
    fn text_editor_accepts_typing_enter_and_escape_without_resetting_or_leaving_fullscreen() {
        let ctx = egui::Context::default();
        let mut app = app();
        let mut frame = eframe::Frame::_new_kittest();
        let key = |key| egui::Event::Key {
            key,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        };
        let mut step = |app: &mut App, events: Vec<egui::Event>| {
            let mut input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(1320., 900.),
                )),
                focused: true,
                events,
                ..Default::default()
            };
            input
                .viewports
                .get_mut(&egui::ViewportId::ROOT)
                .unwrap()
                .fullscreen = Some(true);
            ctx.run_ui(input, |ui| app.ui(ui, &mut frame))
        };
        app.begin_edit(TextTarget::FigureTitle, None);
        app.editor.as_mut().unwrap().value.clear();
        step(&mut app, vec![]).drop_without_applying_deltas();
        step(&mut app, vec![egui::Event::Text("An edited title".into())])
            .drop_without_applying_deltas();
        step(&mut app, vec![key(egui::Key::Enter)]).drop_without_applying_deltas();
        assert_eq!(app.fig.title, "An edited title");
        assert!(app.editor.is_none());
        let views = app.views.clone();
        app.begin_edit(TextTarget::PanelCaption(0), None);
        app.editor.as_mut().unwrap().value = "Discard me".into();
        let output = step(&mut app, vec![key(egui::Key::Escape)]);
        assert!(app.editor.is_none());
        assert!(app.fig.plots[0].caption.is_empty());
        assert_eq!(app.views, views);
        assert!(
            !output.viewport_output[&egui::ViewportId::ROOT]
                .commands
                .iter()
                .any(|c| matches!(c, egui::ViewportCommand::Fullscreen(false)))
        );
        output.drop_without_applying_deltas();
    }
    #[test]
    fn editing_undo_reorder_and_project_reopen_keep_panel_identity() {
        let mut a = app();
        let original = document::Presentation::capture(&a.fig);
        a.begin_edit(TextTarget::PanelCaption(0), None);
        a.editor.as_mut().unwrap().value = "Edited caption α with r^{2}.".into();
        assert!(a.commit_edit());
        assert!(a.modified());
        let edited = document::Presentation::capture(&a.fig);
        a.focus = Some(0);
        a.views[0] = a.views[0].pan(0.2, 0.1);
        let view = a.views[0];
        a.move_panel(0, 2);
        assert_eq!(a.focus, Some(2));
        assert_eq!(a.fig.plots[2].id, 0);
        a.undo(false);
        assert_eq!(a.focus, Some(0));
        assert_eq!(a.views[0], view);
        assert_eq!(document::Presentation::capture(&a.fig), edited);
        a.undo(false);
        assert_eq!(document::Presentation::capture(&a.fig), original);
        a.undo(true);
        assert_eq!(document::Presentation::capture(&a.fig), edited);
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("edited.peteplot.json");
        document::save(&path, a.origin.as_ref().unwrap(), &a.snapshot(), false).unwrap();
        let mut reopened = app();
        reopened.install_project(document::load(&path).unwrap(), path);
        assert_eq!(reopened.snapshot(), a.snapshot());
        assert!(!reopened.modified());
        assert!(!reopened.history.can_undo());
        // No-op commits and cancelled editing do not enter history.
        reopened.begin_edit(TextTarget::FigureTitle, None);
        assert!(reopened.commit_edit());
        assert!(!reopened.history.can_undo());
        reopened.begin_edit(TextTarget::FigureCaption, None);
        reopened.editor = None;
        assert!(!reopened.modified());
    }
    #[test]
    fn text_regions_support_click_editing_at_different_scales_and_in_focused_view() {
        for scale in [1., 2.] {
            for focus in [None, Some(2)] {
                let ctx = egui::Context::default();
                let mut a = app();
                a.focus = focus;
                canvas_frame_scaled(&mut a, &ctx, vec![], scale);
                canvas_frame_scaled(&mut a, &ctx, vec![], scale);
                let region = a
                    .text_regions
                    .iter()
                    .find(|r| r.target == TextTarget::FigureTitle)
                    .unwrap();
                let pos = egui::pos2((region.rect[0] + 10.) as f32, (region.rect[1] + 10.) as f32);
                canvas_frame_scaled(
                    &mut a,
                    &ctx,
                    vec![
                        egui::Event::PointerMoved(pos),
                        egui::Event::PointerButton {
                            pos,
                            button: egui::PointerButton::Primary,
                            pressed: true,
                            modifiers: Default::default(),
                        },
                    ],
                    scale,
                );
                canvas_frame_scaled(
                    &mut a,
                    &ctx,
                    vec![egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed: false,
                        modifiers: Default::default(),
                    }],
                    scale,
                );
                assert_eq!(
                    a.editor.as_ref().unwrap().target,
                    focus
                        .map(TextTarget::PanelTitle)
                        .unwrap_or(TextTarget::FigureTitle)
                );
                assert!(a.drag_panel.is_none());
            }
        }
    }
    #[test]
    fn oversized_caption_does_not_commit_or_lose_the_existing_figure() {
        let mut a = app();
        let before = document::Presentation::capture(&a.fig);
        a.begin_edit(TextTarget::PanelCaption(0), None);
        a.editor.as_mut().unwrap().value = "Long caption with many lines. ".repeat(200);
        assert!(!a.commit_edit());
        assert!(a.editor.is_some());
        assert_eq!(document::Presentation::capture(&a.fig), before);
        assert!(!a.history.can_undo());
        assert!(a.status.contains("would not fit"));
    }
    #[test]
    fn escape_exits_fullscreen_without_closing_the_window() {
        let ctx = egui::Context::default();
        let mut input = egui::RawInput::default();
        input
            .viewports
            .get_mut(&egui::ViewportId::ROOT)
            .unwrap()
            .fullscreen = Some(true);
        input.events.push(egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        let out = ctx.run_ui(input, |ui| exit_fullscreen_on_escape(ui.ctx()));
        let commands = &out.viewport_output[&egui::ViewportId::ROOT].commands;
        assert!(
            commands
                .iter()
                .any(|c| matches!(c, egui::ViewportCommand::Fullscreen(false)))
        );
        assert!(
            !commands
                .iter()
                .any(|c| matches!(c, egui::ViewportCommand::Close))
        );
        out.drop_without_applying_deltas();
    }
    fn canvas_frame(app: &mut App, ctx: &egui::Context, events: Vec<egui::Event>) {
        canvas_frame_scaled(app, ctx, events, 1.);
    }
    fn canvas_frame_scaled(
        app: &mut App,
        ctx: &egui::Context,
        events: Vec<egui::Event>,
        scale: f32,
    ) {
        let mut input = egui::RawInput {
            max_texture_side: Some(4096),
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1200., 960.),
            )),
            events,
            focused: true,
            ..Default::default()
        };
        input
            .viewports
            .get_mut(&egui::ViewportId::ROOT)
            .unwrap()
            .native_pixels_per_point = Some(scale);
        ctx.run_ui(input, |ui| app.canvas(ui))
            .drop_without_applying_deltas();
    }
    #[test]
    fn egui_pointer_events_drive_hover_zoom_and_pan() {
        let ctx = egui::Context::default();
        let mut app = app();
        canvas_frame(&mut app, &ctx, vec![]);
        canvas_frame(&mut app, &ctx, vec![]);
        assert!(!app.error, "{}", app.status);
        assert_eq!(app.areas.len(), 4);
        let a = app.areas[0];
        let pos = egui::pos2((a.x + a.w / 2.) as f32, (a.y + a.h / 2.) as f32);
        canvas_frame(&mut app, &ctx, vec![egui::Event::PointerMoved(pos)]);
        assert!(app.hover.contains("sample"), "{}", app.hover);
        let original = app.views.clone();
        canvas_frame(
            &mut app,
            &ctx,
            vec![egui::Event::MouseWheel {
                unit: egui::MouseWheelUnit::Point,
                delta: egui::vec2(0., 100.),
                phase: egui::TouchPhase::Move,
                modifiers: Default::default(),
            }],
        );
        assert_eq!(
            app.views, original,
            "plain scrolling must not change the view"
        );
        let shift = egui::Modifiers::SHIFT;
        canvas_frame(
            &mut app,
            &ctx,
            vec![
                egui::Event::ModifiersChanged(shift),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: shift,
                },
            ],
        );
        let end = pos + egui::vec2(90., 70.);
        canvas_frame(&mut app, &ctx, vec![egui::Event::PointerMoved(end)]);
        assert_eq!(
            app.views, original,
            "zoom commits on release, not during selection"
        );
        canvas_frame(
            &mut app,
            &ctx,
            vec![
                egui::Event::PointerButton {
                    pos: end,
                    button: egui::PointerButton::Primary,
                    pressed: false,
                    modifiers: shift,
                },
                egui::Event::ModifiersChanged(egui::Modifiers::NONE),
            ],
        );
        assert!(app.views[0].x[1] - app.views[0].x[0] < original[0].x[1] - original[0].x[0]);
        assert!(app.views[0].y[1] - app.views[0].y[0] < original[0].y[1] - original[0].y[0]);
        assert_eq!(
            app.views[1], original[1],
            "zoom must affect only the selected panel"
        );
        canvas_frame(&mut app, &ctx, vec![egui::Event::PointerMoved(pos)]);
        let zoomed = app.views[0];
        canvas_frame(
            &mut app,
            &ctx,
            vec![egui::Event::PointerButton {
                pos,
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: Default::default(),
            }],
        );
        canvas_frame(
            &mut app,
            &ctx,
            vec![egui::Event::PointerMoved(pos + egui::vec2(35., 0.))],
        );
        assert!(
            app.views[0].x[0] < zoomed.x[0],
            "drag should translate the x bounds"
        );
        // Releasing a pan over a title must end the gesture, not open an editor.
        let over_title = egui::pos2(50., 50.);
        canvas_frame(&mut app, &ctx, vec![egui::Event::PointerMoved(over_title)]);
        canvas_frame(
            &mut app,
            &ctx,
            vec![egui::Event::PointerButton {
                pos: over_title,
                button: egui::PointerButton::Primary,
                pressed: false,
                modifiers: Default::default(),
            }],
        );
        assert!(app.drag_panel.is_none());
        assert!(app.editor.is_none());
        app.reset();
        assert_eq!(app.views, original);
    }
    #[test]
    fn save_uses_current_bounds_and_requires_confirmation() {
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("view.svg");
        let mut app = app();
        app.save_path = file.to_str().unwrap().into();
        let original = app.views[0];
        app.views[0] = original.zoom(60., 0.4, 0.5);
        app.save_to_path(file.clone());
        assert!(!app.error, "{}", app.status);
        let zoomed = std::fs::read(&file).unwrap();
        app.full_plot = true;
        app.save_to_path(file.clone());
        let pending = app
            .pending
            .take()
            .expect("existing file must require confirmation");
        assert_eq!(std::fs::read(&file).unwrap(), zoomed);
        assert_eq!(pending.views[0], original);
        app.save(pending, true);
        assert!(!app.error);
        assert_ne!(std::fs::read(&file).unwrap(), zoomed);
        app.full_plot = false;
        app.save_to_path(file.clone());
        let pending = app.pending.take().unwrap();
        assert_eq!(pending.views[0], app.views[0]);
        // Dropping the pending request is cancellation: no write has occurred.
        let before = std::fs::read(&file).unwrap();
        drop(pending);
        assert_eq!(std::fs::read(&file).unwrap(), before);
    }
    #[test]
    fn arrangement_keeps_views_and_selection_with_the_graph() {
        let mut app = app();
        app.focus = Some(0);
        app.views[0] = app.views[0].zoom(60., 0.4, 0.5);
        let title = app.fig.plots[0].title.clone();
        let bounds = app.views[0];
        app.move_panel(0, 1);
        assert_eq!(app.focus, Some(1));
        assert_eq!(app.fig.plots[1].title, title);
        assert_eq!(app.views[1], bounds);
        let (single, views) = app.selected_figure(Some(1));
        assert_eq!(single.plots.len(), 1);
        assert_eq!(single.title, title);
        assert_eq!(views, vec![bounds]);
        let scene = app.renderer.scene(&single, &views, 900, 600).unwrap();
        assert_eq!(scene.areas.len(), 1);
        assert_eq!(scene.areas[0].bounds, bounds);
        let (all, views) = app.selected_figure(None);
        assert_eq!(all.plots.len(), 4);
        assert_eq!(views[1], bounds);
    }
    #[test]
    fn thumbnail_click_and_drag_select_and_move_graphs() {
        let ctx = egui::Context::default();
        let mut app = app();
        let draw = |app: &mut App, events| {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(200., 1500.),
                )),
                events,
                ..Default::default()
            };
            ctx.run_ui(input, |ui| app.sidebar(ui))
                .drop_without_applying_deltas();
        };
        draw(&mut app, vec![]);
        draw(&mut app, vec![]);
        let first = ctx
            .read_response(egui::Id::new(("graph", 0usize)))
            .unwrap()
            .rect
            .center();
        let third = ctx
            .read_response(egui::Id::new(("graph", 2usize)))
            .unwrap()
            .rect
            .center();
        let button = |pos, pressed| egui::Event::PointerButton {
            pos,
            pressed,
            button: egui::PointerButton::Primary,
            modifiers: Default::default(),
        };
        draw(
            &mut app,
            vec![egui::Event::PointerMoved(third), button(third, true)],
        );
        draw(&mut app, vec![button(third, false)]);
        assert_eq!(app.focus, Some(2));
        let original: Vec<_> = app.fig.plots.iter().map(|p| p.title.clone()).collect();
        draw(&mut app, vec![button(third, true)]);
        draw(
            &mut app,
            vec![egui::Event::PointerMoved(third + egui::vec2(0., -20.))],
        );
        draw(&mut app, vec![egui::Event::PointerMoved(first)]);
        draw(&mut app, vec![button(first, false)]);
        assert_eq!(app.focus, Some(0));
        let reordered: Vec<_> = app.fig.plots.iter().map(|p| p.title.clone()).collect();
        assert_eq!(
            reordered,
            vec![
                original[2].clone(),
                original[0].clone(),
                original[1].clone(),
                original[3].clone()
            ]
        );
    }
    #[test]
    fn right_click_selects_a_panel_without_panning() {
        let ctx = egui::Context::default();
        let mut app = app();
        canvas_frame(&mut app, &ctx, vec![]);
        canvas_frame(&mut app, &ctx, vec![]);
        let a = app.areas[2];
        let pos = egui::pos2((a.x + a.w / 2.) as f32, (a.y + a.h / 2.) as f32);
        let original = app.views.clone();
        canvas_frame(
            &mut app,
            &ctx,
            vec![
                egui::Event::PointerMoved(pos),
                egui::Event::PointerButton {
                    pos,
                    button: egui::PointerButton::Secondary,
                    pressed: true,
                    modifiers: Default::default(),
                },
            ],
        );
        canvas_frame(
            &mut app,
            &ctx,
            vec![egui::Event::PointerButton {
                pos,
                button: egui::PointerButton::Secondary,
                pressed: false,
                modifiers: Default::default(),
            }],
        );
        assert_eq!(app.context_panel, Some(2));
        assert_eq!(app.views, original);
    }
    #[test]
    fn full_ui_lays_out_at_supported_window_sizes() {
        let ctx = egui::Context::default();
        let mut app = app();
        let mut frame = eframe::Frame::_new_kittest();
        for size in [[1320., 900.], [850., 600.]] {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
                ..Default::default()
            };
            let out = ctx.run_ui(input, |ui| app.ui(ui, &mut frame));
            assert!(!out.shapes.is_empty());
            out.drop_without_applying_deltas();
            assert!(!app.error, "{}", app.status);
            assert!(app.texture.is_some());
        }
    }
    #[test]
    fn qc_annotations_survive_selection_reorder_and_current_view_export() {
        let args = [
            "-infiles",
            "1D: 0 1 2 3 | 3 2 1 0",
            "-shade",
            "0.5:1.5",
            "-shade-label",
            "Task A",
            "-censor_trs",
            "1",
            "-boxplot_on",
            "-yfiles_pm",
            "1D: 4@0.2 | 4@0.3",
        ];
        let o = crate::cli::parse(&args.map(String::from)).unwrap();
        let f = crate::model::from_options(&o).unwrap();
        let mut app = App::new(f, o, Renderer::new());
        app.views[0] = app.views[0].zoom(1., 1., 0.5);
        let view = app.views[0];
        app.move_panel(0, 1);
        let (f, v) = app.selected_figure(Some(1));
        assert_eq!(v, [view]);
        assert_eq!(f.plots[0].qc.censored, [false, true, false, false]);
        assert_eq!(
            f.plots[0].qc.annotations[0].label.as_deref(),
            Some("Task A")
        );
        assert_eq!(f.plots[0].qc.boxes[1].n, 3);
        assert!(f.plots[0].hover(1., 1.).contains("CENSORED"));
        assert!(app.renderer.scene(&f, &v, 900, 600).is_ok());
    }
    #[test]
    fn attached_graph_add_undo_duplicate_and_replay_keep_scientific_roles() {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/attachments");
        let table = crate::import::inspect(&base.join("signals.tsv"), || false).unwrap();
        let mut m = crate::import::Mapping::for_kind(&table, crate::import::Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.uncertainty = "bounds".into();
        m.attachments.bands.insert(1, [Some(2), Some(3)]);
        m.attachments.group = Some(5);
        m.attachments.censor_column = Some(6);
        m.attachments.xmat = base.join("X.mat.1D").display().to_string();
        m.attachments.aligned = true;
        let preview = crate::import::prepare(&table, &m, || false).unwrap();
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        app.add_import(preview).unwrap();
        app.undo(false);
        assert!(app.fig.plots.is_empty());
        app.undo(true);
        app.duplicate_graph(0).unwrap();
        let pixels = |app: &App| {
            app.renderer
                .raster(
                    &app.renderer
                        .scene(&app.fig, &app.views, app.export_width, app.export_height)
                        .unwrap(),
                )
                .unwrap()
                .data()
                .to_vec()
        };
        let before = pixels(&app);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("attachments.json");
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        app.install_project(document::load(&path).unwrap(), path);
        assert_eq!(pixels(&app), before);
        assert_eq!(app.fig.plots[0].qc.censored, app.fig.plots[1].qc.censored);
        assert_eq!(app.fig.plots[0].qc.annotations.len(), 3);
    }
}
