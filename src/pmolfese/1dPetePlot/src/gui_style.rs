//! Full-window styling workspace with a rendered draft and always-visible commit controls.
use super::*;
use crate::{
    model::Content,
    style::{PALETTES, PlotStyle},
};
use anyhow::Context;

pub(super) struct Editor {
    index: usize,
    draft: PlotStyle,
    shared: PlotStyle,
    all: bool,
    preview: Option<egui::TextureHandle>,
    preview_key: Option<(PlotStyle, PlotStyle, bool)>,
    error: String,
}
impl App {
    pub(super) fn begin_style(&mut self, index: usize) {
        if self.editor.is_some() || index >= self.fig.plots.len() {
            return;
        }
        self.styling = Some(Editor {
            index,
            draft: self.fig.plots[index].style.clone(),
            shared: PlotStyle::default(),
            all: false,
            preview: None,
            preview_key: None,
            error: String::new(),
        });
        self.drag_panel = None;
        self.zoom_start = None;
    }
    fn styled_figure(&self, editor: &Editor) -> Result<Figure> {
        let mut figure = self.fig.clone();
        if editor.all {
            for plot in &mut figure.plots {
                if let Some(width) = editor.shared.line_width
                    && matches!(plot.content, Content::Lines(_) | Content::Intervals { .. })
                {
                    plot.style.line_width = Some(width);
                    for s in &mut plot.style.series {
                        s.line_width = None;
                    }
                }
                if let Some(radius) = editor.shared.marker_radius
                    && matches!(
                        plot.content,
                        Content::Lines(_) | Content::Scatter(_) | Content::Intervals { .. }
                    )
                {
                    plot.style.marker_radius = Some(radius);
                    for s in &mut plot.style.series {
                        s.marker_radius = None;
                    }
                }
                if editor.shared.legend.is_some()
                    && matches!(
                        plot.content,
                        Content::Lines(_) | Content::Scatter(_) | Content::Histogram { .. }
                    )
                {
                    plot.style.legend = editor.shared.legend;
                }
                if matches!(plot.content, Content::Heat { .. }) && editor.shared.palette.is_some() {
                    plot.style.palette.clone_from(&editor.shared.palette);
                }
            }
        } else {
            figure.plots[editor.index].style = editor.draft.clone();
        }
        for plot in &mut figure.plots {
            while plot
                .style
                .series
                .last()
                .is_some_and(|s| *s == Default::default())
            {
                plot.style.series.pop();
            }
            plot.style.validate(plot)?;
        }
        Ok(figure)
    }
    pub(super) fn commit_style(&mut self) -> bool {
        let Some(editor) = &self.styling else {
            return true;
        };
        let result = self.styled_figure(editor).and_then(|fig| {
            self.renderer
                .scene(&fig, &self.views, self.export_width, self.export_height)?;
            Ok(fig)
        });
        match result {
            Ok(figure) => {
                let before = document::Presentation::capture(&self.fig);
                self.fig = figure;
                self.styling = None;
                self.edited(before);
                self.status =
                    "Style applied · Undo available · included in images and saved projects".into();
                self.error = false;
                true
            }
            Err(e) => {
                self.styling.as_mut().unwrap().error = format!(
                    "Style was not applied: {e:#}. Shorten legend labels or cancel and increase image height."
                );
                self.status = "Style could not be applied; return to the style workspace".into();
                self.error = true;
                false
            }
        }
    }
    pub(super) fn style_workspace(&mut self, ui: &mut egui::Ui) {
        let mut editor = self.styling.take().expect("open styling workspace");
        let mut cancel = ui.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape));
        let mut apply = false;
        let plot = &self.fig.plots[editor.index];
        egui::Panel::bottom("style-commit")
            .frame(egui::Frame::new().fill(Color32::WHITE).inner_margin(16.))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    apply |= ui
                        .add_sized(
                            [150., 36.],
                            egui::Button::new(
                                RichText::new(if editor.all {
                                    "Apply to all graphs"
                                } else {
                                    "Apply to this graph"
                                })
                                .color(Color32::WHITE),
                            )
                            .fill(Color32::from_rgb(20, 125, 146)),
                        )
                        .clicked();
                    cancel |= ui
                        .add_sized([90., 36.], egui::Button::new("Cancel"))
                        .clicked();
                    ui.label("One Undo restores the previous styling. Escape cancels.");
                });
                if !editor.error.is_empty() {
                    ui.colored_label(Color32::from_rgb(175, 50, 50), &editor.error);
                }
            });
        egui::Panel::left("style-controls").default_size(335.).size_range(300.0..=500.0).resizable(true).frame(egui::Frame::new().fill(Color32::WHITE).inner_margin(20.)).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("Style your figure");
                ui.label(format!("Graph {} · {}", editor.index + 1, plot.title));
                ui.add_space(8.);
                ui.radio_value(&mut editor.all, false, "This graph");
                ui.radio_value(&mut editor.all, true, "All graphs · shared appearance");
                ui.small(if editor.all { "Checked settings replace existing overrides in all compatible graphs. Colors and legend names stay individual." } else { "Unchecked settings use the original graph defaults." });
                let draft = if editor.all { &mut editor.shared } else { &mut editor.draft };
                ui.separator();
                if editor.all || matches!(plot.content, Content::Lines(_) | Content::Intervals { .. }) {
                    numeric(ui, "Line width", &mut draft.line_width, 2.1, 0.25..=12.);
                }
                if editor.all || matches!(plot.content, Content::Lines(_) | Content::Scatter(_) | Content::Intervals { .. }) {
                    numeric(ui, "Marker radius", &mut draft.marker_radius, 3.4, 0.5..=20.);
                    ui.small("Sizes scale with the image. Lines show markers at isolated observations only.");
                }
                if editor.all || matches!(plot.content, Content::Lines(_) | Content::Scatter(_) | Content::Histogram { .. }) {
                    let mut custom = draft.legend.is_some();
                    if ui.checkbox(&mut custom, "Override series legend").changed() { draft.legend = custom.then_some(plot.legend); }
                    if let Some(show) = &mut draft.legend { ui.checkbox(show, "Show series legend"); }
                    ui.small("Shading keys and reference labels remain visible.");
                }
                if editor.all || matches!(plot.content, Content::Heat { .. }) {
                    let mut custom = draft.palette.is_some();
                    if ui.checkbox(&mut custom, "Override heatmap palette").changed() { draft.palette = custom.then(|| "viridis".into()); }
                    if let Some(palette) = &mut draft.palette {
                        egui::ComboBox::from_id_salt("style-palette").selected_text(palette.clone()).show_ui(ui, |ui| {
                            for name in PALETTES { ui.selectable_value(palette, name.into(), name); }
                        });
                    }
                    ui.small("Color limits and condition ordering are preserved.");
                }
                if !editor.all {
                    ui.separator();
                    let series = plot.style_series();
                    for (index, (label, original)) in series.iter().enumerate() {
                        ui.push_id(index, |ui| {
                            egui::CollapsingHeader::new(format!("Series {} · {}", index + 1, label)).default_open(index == 0).show(ui, |ui| {
                                draft.series.resize_with(series.len(), Default::default);
                                let s = &mut draft.series[index];
                                ui.horizontal(|ui| {
                                    let c = s.color.as_deref().unwrap_or(original);
                                    let mut rgb = [0; 3];
                                    for (i, channel) in rgb.iter_mut().enumerate() { *channel = u8::from_str_radix(&c[1+i*2..3+i*2],16).unwrap_or(0); }
                                    if ui.color_edit_button_srgb(&mut rgb).changed() { s.color = Some(format!("#{:02x}{:02x}{:02x}",rgb[0],rgb[1],rgb[2])); }
                                    ui.label("Color");
                                    if ui.small_button("Original").clicked() { s.color = None; }
                                });
                                if !matches!(plot.content, Content::Intervals { .. }) {
                                    let mut custom = s.label.is_some();
                                    if ui.checkbox(&mut custom, "Custom legend name").changed() { s.label = custom.then(|| label.to_string()); }
                                    if let Some(label) = &mut s.label { ui.add(egui::TextEdit::singleline(label).desired_width(f32::INFINITY)); }
                                }
                                if matches!(plot.content, Content::Lines(_) | Content::Intervals { .. }) {
                                    numeric(ui, "Series line width", &mut s.line_width, draft.line_width.unwrap_or(2.1), 0.25..=12.);
                                }
                                if matches!(plot.content, Content::Lines(_) | Content::Scatter(_) | Content::Intervals { .. }) {
                                    numeric(ui, "Series marker radius", &mut s.marker_radius, draft.marker_radius.unwrap_or(3.4), 0.5..=20.);
                                }
                            });
                        });
                    }
                    ui.add_space(12.);
                    if ui.button("Restore original graph styling").clicked() { *draft = PlotStyle::default(); }
                }
            });
        });
        let key = (editor.draft.clone(), editor.shared.clone(), editor.all);
        if editor.preview_key.as_ref() != Some(&key) {
            editor.preview_key = Some(key);
            editor.error.clear();
            let preview = self.styled_figure(&editor).and_then(|mut fig| {
                let views = if editor.all {
                    self.views.clone()
                } else {
                    let p = fig.plots.remove(editor.index);
                    fig.title = p.title.clone();
                    fig.caption.clear();
                    fig.columns = 1;
                    fig.plots = vec![p];
                    vec![self.views[editor.index]]
                };
                let (w, h) = if editor.all {
                    (
                        720,
                        (self.export_height as f64 * 720. / self.export_width as f64).round()
                            as u32,
                    )
                } else {
                    (720, 600)
                };
                let scene = self.renderer.scene(&fig, &views, w, h)?;
                self.renderer.raster(&scene)
            });
            match preview {
                Ok(pixels) => {
                    editor.preview = Some(ui.ctx().load_texture(
                        "style-preview",
                        egui::ColorImage::from_rgba_premultiplied(
                            [pixels.width() as usize, pixels.height() as usize],
                            pixels.data(),
                        ),
                        egui::TextureOptions::LINEAR,
                    ));
                }
                Err(e) => {
                    editor.preview = None;
                    editor.error = format!("Preview: {e:#}");
                }
            }
        }
        egui::CentralPanel::default()
            .frame(
                egui::Frame::new()
                    .fill(Color32::from_rgb(237, 242, 244))
                    .inner_margin(20.),
            )
            .show(ui, |ui| {
                ui.label(if editor.all {
                    "PREVIEW · ALL GRAPHS"
                } else {
                    "PREVIEW · SELECTED GRAPH"
                });
                if let Some(texture) = &editor.preview {
                    let available = ui.available_size();
                    ui.add(
                        egui::Image::new(texture).fit_to_exact_size(
                            texture.size_vec2()
                                * (available.x / texture.size_vec2().x)
                                    .min(available.y / texture.size_vec2().y)
                                    .max(0.01),
                        ),
                    );
                }
            });
        self.styling = Some(editor);
        if cancel {
            self.styling = None;
        } else if apply {
            self.commit_style();
        }
    }
    pub(super) fn duplicate_graph(&mut self, index: usize) -> Result<()> {
        anyhow::ensure!(
            self.fig.plots.len() + self.dormant.len() < 64,
            "This session retains 64 graphs including Undo history. Save and reopen before adding more."
        );
        anyhow::ensure!(
            self.fig.plots.iter().all(|p| !p.qc.margin_off),
            "Margin-free QC strips cannot share a layout. Reopen without -margin_off to duplicate."
        );
        let original = self
            .fig
            .plots
            .get(index)
            .context("Select a graph to duplicate")?;
        let source = original.id;
        let mut copy = original.clone();
        copy.id = self.next_id;
        copy.title = format!("{} · copy", original.title);
        let mut figure = self.fig.clone();
        let mut views = self.views.clone();
        figure.plots.insert(index + 1, copy);
        views.insert(index + 1, views[index]);
        let minimum = ((150. + figure.plots.len().div_ceil(figure.columns) as f64 * 290.)
            * self.export_width as f64
            / 1200.)
            .ceil() as u32;
        let mut height = self.export_height.max(minimum);
        let max_height = 8192.min(32_000_000 / self.export_width);
        loop {
            anyhow::ensure!(
                height <= max_height,
                "The duplicate exceeds image limits. Increase Columns or reduce export width, then retry."
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
        let before = document::Presentation::capture(&self.fig);
        if let Some(origin) = &mut self.origin {
            origin.duplicate(self.next_id, source);
        }
        self.next_id += 1;
        self.fig = figure;
        self.views = views;
        self.export_height = height;
        self.focus = Some(index + 1);
        self.edit_panel = index + 1;
        self.edited(before);
        self.status = "Graph duplicated · edit the copy independently · Undo available".into();
        self.error = false;
        Ok(())
    }
}
fn numeric(
    ui: &mut egui::Ui,
    label: &str,
    value: &mut Option<f64>,
    default: f64,
    range: std::ops::RangeInclusive<f64>,
) {
    ui.horizontal(|ui| {
        let mut custom = value.is_some();
        if ui.checkbox(&mut custom, label).changed() {
            *value = custom.then_some(default);
        }
        if let Some(value) = value {
            ui.add(egui::DragValue::new(value).range(range).speed(0.1));
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::style::SeriesStyle;
    use eframe::App as _;
    fn app() -> App {
        App::new(crate::model::demo(), Options::default(), Renderer::new())
    }
    fn pixels(app: &App) -> Vec<u8> {
        app.renderer
            .raster(
                &app.renderer
                    .scene(&app.fig, &app.views, app.export_width, app.export_height)
                    .unwrap(),
            )
            .unwrap()
            .data()
            .to_vec()
    }
    #[test]
    fn selected_and_shared_styling_preserve_data_and_are_atomic_undo_actions() {
        let mut app = app();
        let before = app.snapshot();
        let original_pixels = pixels(&app);
        let source = format!("{:?}", app.fig.plots[0].content);
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.series = vec![SeriesStyle {
            color: Some("#bc477b".into()),
            label: Some("Custom response".into()),
            line_width: Some(5.),
            marker_radius: None,
        }];
        assert!(app.commit_style());
        assert_ne!(pixels(&app), original_pixels);
        assert_eq!(format!("{:?}", app.fig.plots[0].content), source);
        assert_eq!(app.views, before.views);
        for i in 1..app.fig.plots.len() {
            assert_eq!(app.fig.plots[i].style, PlotStyle::default());
        }
        let styled = app.snapshot();
        app.undo(false);
        assert_eq!(app.snapshot(), before);
        app.undo(true);
        assert_eq!(app.snapshot(), styled);
        app.begin_style(0);
        let editor = app.styling.as_mut().unwrap();
        editor.all = true;
        editor.shared = PlotStyle {
            line_width: Some(4.),
            marker_radius: Some(6.),
            legend: Some(false),
            palette: Some("Reds_and_Blues_Inv".into()),
            ..Default::default()
        };
        assert!(app.commit_style());
        assert_eq!(app.fig.plots[0].style.width(0, 2.1), 4.);
        assert_eq!(
            app.fig.plots[0].style.label(0, "original"),
            "Custom response"
        );
        assert_eq!(app.fig.plots[0].style.color(0, "original"), "#bc477b");
        assert_eq!(app.fig.plots[1].style.line_width, None); // histogram has no line stroke
        assert_eq!(app.fig.plots[2].style.radius(0, 3.4), 6.);
        assert_eq!(app.fig.plots[3].style.cmap("viridis"), "Reds_and_Blues_Inv");
        app.undo(false);
        assert_eq!(app.snapshot(), styled); // all changes undone together
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.series[0].color = Some("<invalid>".into());
        assert!(!app.commit_style());
        assert_eq!(app.snapshot(), styled);
        assert!(app.styling.is_some());
    }
    #[test]
    fn duplicate_copy_reorder_branch_and_replay_preserve_independent_styles_and_views() {
        let mut app = app();
        app.views[0] = app.views[0].pan(0.5, 0.2);
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.line_width = Some(3.);
        assert!(app.commit_style());
        app.duplicate_graph(0).unwrap();
        let first_copy_id = app.fig.plots[1].id;
        assert_eq!(app.views[0], app.views[1]);
        assert_eq!(app.fig.plots[0].style, app.fig.plots[1].style);
        assert_eq!(
            format!("{:?}", app.fig.plots[0].content),
            format!("{:?}", app.fig.plots[1].content)
        );
        app.begin_style(1);
        app.styling.as_mut().unwrap().draft.series = vec![SeriesStyle {
            color: Some("#e18a47".into()),
            label: Some("Independent copy".into()),
            ..Default::default()
        }];
        assert!(app.commit_style());
        assert_ne!(app.fig.plots[0].style, app.fig.plots[1].style);
        app.duplicate_graph(1).unwrap();
        let second_copy_id = app.fig.plots[2].id;
        app.move_panel(2, 0);
        assert_eq!(app.fig.plots[0].id, second_copy_id);
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("styled.json");
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), false).unwrap();
        let loaded = document::load(&path).unwrap();
        let expected = pixels(&app);
        app.install_project(loaded, path.clone());
        assert_eq!(pixels(&app), expected);
        assert_eq!(app.fig.plots[2].id, first_copy_id);
        app.duplicate_graph(0).unwrap();
        app.undo(false);
        assert_eq!(app.fig.plots.len(), 6);
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), true).unwrap();
        assert_eq!(document::load(&path).unwrap().figure.plots.len(), 6);
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.line_width = Some(7.);
        assert!(app.commit_style());
        assert!(!app.history.can_redo());
        assert!(app.dormant.is_empty());
        document::save(&path, app.origin.as_ref().unwrap(), &app.snapshot(), true).unwrap();
        assert_eq!(document::load(&path).unwrap().figure.plots.len(), 6);
    }
    #[test]
    fn imported_ci_and_rdm_duplicates_replay_without_changing_scientific_content() {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let mut app = App::new(Figure::empty(), Options::default(), Renderer::new());
        for (file, kind) in [
            (
                "examples/rsa/ordinary.rsa.1D",
                crate::import::Kind::Intervals,
            ),
            ("examples/rdm/model.1D", crate::import::Kind::Rdm),
        ] {
            let table = crate::import::inspect(&base.join(file), || false).unwrap();
            let mut mapping = crate::import::Mapping::for_kind(&table, kind);
            if kind == crate::import::Kind::Rdm {
                mapping.rdm.labels = base
                    .join("examples/rdm/conditions.txt")
                    .to_str()
                    .unwrap()
                    .into();
            }
            let prepared = crate::import::prepare(&table, &mapping, || false).unwrap();
            app.add_import(prepared).unwrap();
            let i = app.fig.plots.len() - 1;
            app.duplicate_graph(i).unwrap();
            app.begin_style(i + 1);
            let draft = &mut app.styling.as_mut().unwrap().draft;
            if kind == crate::import::Kind::Intervals {
                draft.series = vec![SeriesStyle {
                    color: Some("#7064b5".into()),
                    ..Default::default()
                }];
                draft.line_width = Some(4.);
                draft.marker_radius = Some(6.);
            } else {
                draft.palette = Some("Reds_and_Blues".into());
            }
            assert!(app.commit_style());
            assert_eq!(
                format!("{:?}", app.fig.plots[i].content),
                format!("{:?}", app.fig.plots[i + 1].content)
            );
        }
        let path = tempfile::tempdir().unwrap();
        let project = path.path().join("mixed.json");
        let expected = pixels(&app);
        document::save(
            &project,
            app.origin.as_ref().unwrap(),
            &app.snapshot(),
            false,
        )
        .unwrap();
        app.install_project(document::load(&project).unwrap(), project);
        assert_eq!(pixels(&app), expected);
    }
    #[test]
    fn full_window_style_footer_applies_without_scrolling_and_escape_cancels_in_fullscreen() {
        let mut app = app();
        let ctx = egui::Context::default();
        let mut frame = eframe::Frame::_new_kittest();
        let draw = |app: &mut App, events| {
            let mut input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    egui::vec2(850., 500.),
                )),
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
        let mut draw = draw;
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.line_width = Some(5.);
        draw(&mut app, vec![]).drop_without_applying_deltas();
        draw(&mut app, vec![]).drop_without_applying_deltas();
        // Fixed footer: left action stays at the bottom even with a scrollable inspector.
        let pos = egui::pos2(80., 465.);
        for pressed in [true, false] {
            draw(
                &mut app,
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    },
                ],
            )
            .drop_without_applying_deltas();
        }
        assert!(
            app.styling.is_none(),
            "Apply should be visible at 850x500 without scrolling"
        );
        assert_eq!(app.fig.plots[0].style.line_width, Some(5.));
        let before = app.snapshot();
        app.begin_style(0);
        app.styling.as_mut().unwrap().draft.line_width = Some(9.);
        let output = draw(
            &mut app,
            vec![egui::Event::Key {
                key: egui::Key::Escape,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: Default::default(),
            }],
        );
        assert!(app.styling.is_none());
        assert_eq!(app.snapshot(), before);
        assert!(
            !output.viewport_output[&egui::ViewportId::ROOT]
                .commands
                .iter()
                .any(|c| matches!(c, egui::ViewportCommand::Fullscreen(false)))
        );
        output.drop_without_applying_deltas();
    }
}
