//! Native-window import workflow. Temporary candidates never mutate the document.
use crate::import::{Answer, Kind, Mapping, Prepared, Table, Work, Worker};
use eframe::egui::{self, Color32, RichText};
use std::{collections::VecDeque, path::PathBuf, sync::Arc};

#[derive(Default)]
pub struct ImportUi {
    worker: Option<Worker>,
    queue: VecDeque<PathBuf>,
    path: Option<PathBuf>,
    table: Option<Arc<Table>>,
    mapping: Option<Mapping>,
    preview: Option<Box<Prepared>>,
    texture: Option<egui::TextureHandle>,
    busy: bool,
    pub error: Option<String>,
    #[cfg(test)]
    add_rect: Option<egui::Rect>,
    #[cfg(test)]
    pub(crate) page_rect: Option<egui::Rect>,
}
impl ImportUi {
    pub fn visible(&self) -> bool {
        self.path.is_some()
    }
    pub fn enqueue(&mut self, paths: Vec<PathBuf>, ctx: &egui::Context) {
        if paths.len() + self.queue.len() + usize::from(self.visible()) > 16 {
            self.error = Some("Import up to 16 files at a time.".into());
            return;
        }
        if self.worker.is_none() {
            let ctx = ctx.clone();
            self.worker = Some(Worker::new(move || ctx.request_repaint()));
        }
        self.queue.extend(paths);
        if self.path.is_none() {
            self.advance();
        }
    }
    pub fn cancel(&mut self) {
        if let Some(w) = &mut self.worker {
            w.cancel();
        }
        self.queue.clear();
        self.path = None;
        self.table = None;
        self.mapping = None;
        self.preview = None;
        self.texture = None;
        self.busy = false;
        self.error = None;
    }
    pub fn advance(&mut self) {
        self.path = self.queue.pop_front();
        self.retry();
    }
    fn retry(&mut self) {
        self.table = None;
        self.mapping = None;
        self.preview = None;
        self.texture = None;
        self.error = None;
        self.busy = self.path.is_some();
        if let (Some(worker), Some(path)) = (&mut self.worker, &self.path) {
            worker.start(Work::Inspect(path.clone()));
        }
    }
    fn refresh(&mut self) {
        self.preview = None;
        self.texture = None;
        self.error = None;
        if let (Some(worker), Some(table), Some(mapping)) =
            (&mut self.worker, &self.table, &self.mapping)
        {
            self.busy = true;
            worker.start(Work::Preview(table.clone(), mapping.clone()));
        }
    }
    pub fn reject(&mut self, message: String) {
        self.error = Some(message);
        self.busy = false;
    }
    pub fn show(&mut self, ui: &mut egui::Ui) -> Option<Box<Prepared>> {
        let ctx = ui.ctx().clone();
        if !self.visible() {
            return None;
        }
        if self.busy {
            ctx.request_repaint_after(std::time::Duration::from_millis(50));
        }
        if let Some(answer) = self.worker.as_mut().and_then(Worker::poll) {
            self.busy = false;
            match answer {
                Ok(Answer::Table(table)) => {
                    self.mapping = Some(Mapping::for_kind(
                        &table,
                        if table.rsa.is_some() {
                            Kind::Histogram
                        } else {
                            Kind::Line
                        },
                    ));
                    self.table = Some(table);
                    self.refresh();
                }
                Ok(Answer::Preview(preview)) => {
                    self.texture = Some(ctx.load_texture(
                        "import-preview",
                        egui::ColorImage::from_rgba_premultiplied(preview.size, &preview.rgba),
                        egui::TextureOptions::LINEAR,
                    ));
                    if self.mapping.as_ref() == Some(&preview.mapping) {
                        self.preview = Some(preview);
                    }
                }
                Ok(Answer::Accepted(prepared)) => return Some(prepared),
                Err(message) => self.error = Some(message),
            }
        }
        let mut close = false;
        let mut skip = false;
        let mut retry = false;
        let mut changed = false;
        let mut add = false;
        let path = self.path.as_ref().unwrap();
        let filename = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        let page = egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(Color32::from_rgb(243, 246, 247)).inner_margin(18.))
            .show(ui, |ui| {
                // Reserve the footer before laying out content: actions must never scroll away.
                let body_height = (ui.available_height() - 48.).max(100.);
                egui::ScrollArea::vertical().id_salt("import-body").max_height(body_height).auto_shrink([false, false]).show(ui, |ui| {
                ui.label(RichText::new("IMPORT DATA").size(11.).color(Color32::from_rgb(20,125,146)));
                ui.heading(&filename);
                ui.label("Preview your graph, then add it to the figure. Cancel returns to your figure.");
                ui.small(path.display().to_string());
                if !self.queue.is_empty() { ui.label(format!("{} more file(s) queued · each gets its own preview", self.queue.len())); }
                if let (Some(table), Some(mapping)) = (&self.table, &mut self.mapping) {
                    ui.label(format!("{} rows · {} columns · {} {}", table.rows, table.columns.len(), table.numeric().len(), if table.rsa.is_some() { "RSA effects" } else { "numeric" }));
                    if table.rsa.is_some() {
                        ui.label("3dRSA results · choose effects of the same kind (r, b, or pr). Statistics are retained as metadata.");
                        ui.horizontal(|ui| { ui.label("ROI IDs (blank = all)"); changed |= ui.text_edit_singleline(&mut mapping.roi_filter).on_hover_text("Comma-separated ROI IDs, e.g. 3,12,27. The chosen order is retained.").changed(); });
                    }
                    ui.horizontal_wrapped(|ui| {
                        for kind in Kind::ALL {
                            let enabled = match kind { Kind::Line if table.rsa.is_some() => false, Kind::Scatter => table.rsa.as_ref().map(|rsa| ["r", "b", "pr"].iter().any(|kind| rsa.effects.keys().filter(|(_, k)| k == kind).count() >= 2)).unwrap_or_else(|| table.numeric().len() >= 2), Kind::Heat => table.rsa.is_some() || table.square(), _ => !table.numeric().is_empty() };
                            if ui.add_enabled(enabled, egui::Button::new(kind.name()).selected(mapping.kind == kind)).on_hover_text(kind.explanation()).clicked() {
                                let missing = mapping.missing.clone(); let rois = mapping.roi_filter.clone(); *mapping = Mapping::for_kind(table, kind); mapping.missing = missing; mapping.roi_filter = rois; changed = true;
                            }
                        }
                    });
                    let explanation = if table.rsa.is_some() { match mapping.kind {
                        Kind::Histogram => "One model effect across selected ROIs. Add other model distributions separately.",
                        Kind::Scatter => "Two effects of the same kind, paired by ROI identity; hover for ROI labels.",
                        Kind::Heat => "ROI-by-model effects (up to 40 selected ROIs), not an RDM. Gray cells are missing.",
                        _ => "Choose an RSA plot type.",
                    }} else { mapping.kind.explanation() };
                    ui.label(RichText::new(explanation).color(Color32::from_rgb(20,125,146)));
                    ui.separator();
                    ui.horizontal_top(|ui| {
                        ui.vertical(|ui| {
                            ui.set_width(235.);
                            if matches!(mapping.kind, Kind::Line | Kind::Scatter) {
                                ui.label("X coordinates");
                                egui::ComboBox::from_id_salt("import-x").selected_text(mapping.x.map(|i| table.columns[i].name.as_str()).unwrap_or("Row index (0-based)"))
                                    .width(225.).show_ui(ui, |ui| {
                                        if mapping.kind == Kind::Line { changed |= ui.selectable_value(&mut mapping.x, None, "Row index (0-based)").changed(); }
                                        for i in table.numeric() { changed |= ui.selectable_value(&mut mapping.x, Some(i), &table.columns[i].name).changed(); }
                                    });
                            }
                            ui.label(if table.rsa.is_some() { "Y / model effect columns" } else if mapping.kind == Kind::Heat { "Matrix columns (row order retained)" } else { "Y / value columns" });
                            egui::ScrollArea::vertical().id_salt("import-columns").max_height(190.).show(ui, |ui| {
                                for (i, column) in table.columns.iter().enumerate() {
                                    let mut selected = mapping.y.contains(&i);
                                    let label = if !column.role.is_empty() { format!("{} · {}", column.name, column.role) } else { format!("{}{}", column.name, if column.invalid { " · infinity" } else if !column.numeric { " · text" } else { "" }) };
                                    let response = ui.add_enabled(column.numeric && !column.invalid, egui::Checkbox::new(&mut selected, label));
                                    if response.changed() {
                                        if selected { mapping.y.push(i); mapping.y.sort_unstable(); } else { mapping.y.retain(|&c| c != i); } changed = true;
                                    }
                                    if column.missing > 0 { ui.small(format!("{} missing", column.missing)); }
                                }
                            });
                            ui.separator(); ui.label("Missing observations");
                            egui::ComboBox::from_id_salt("import-missing").selected_text(match mapping.missing.as_str() { "gap" => "Allow gaps / omit from counts", "mean" => "Impute column mean", _ => "Reject missing values" })
                                .width(225.).show_ui(ui, |ui| {
                                    for (value,label) in [("error","Reject missing values"),("gap","Allow gaps / omit from counts"),("mean","Impute column mean")] {
                                        if table.rsa.is_some() && value == "mean" { continue; }
                                        changed |= ui.selectable_value(&mut mapping.missing, value.into(), label).changed();
                                    }
                                });
                            ui.small(if table.rsa.is_some() { "Missing scatter pairs are omitted together with gap policy. No significance filtering or imputation." } else { "Text columns stay in the table; they are not converted to numbers. X must be finite." });
                            if mapping.missing == "mean" { ui.colored_label(Color32::from_rgb(165,95,20), "Mean imputation changes selected missing values and is recorded in the project."); }
                        });
                        ui.separator();
                        ui.vertical(|ui| {
                            ui.set_min_width(260.);
                            if let Some(texture) = &self.texture {
                                let aspect = texture.size()[0] as f32 / texture.size()[1] as f32;
                                let w = ui.available_width().min((body_height - 190.).max(180.) * aspect); ui.add(egui::Image::new((texture.id(), egui::vec2(w, w / aspect))));
                                ui.small("Preview uses all selected observations. No downsampling.");
                            } else { ui.add_space(50.); ui.label(if self.busy { "Preparing preview…" } else { "Choose compatible columns to preview this graph." }); }
                        });
                    });
                    ui.collapsing("Inspect source rows (first 8)", |ui| {
                        egui::ScrollArea::both().id_salt("source-rows").max_height(165.).show(ui, |ui| {
                            egui::Grid::new("import-grid").striped(true).show(ui, |ui| {
                                ui.strong("Row"); for column in &table.columns { ui.strong(&column.name); } ui.end_row();
                                for (i,row) in table.preview.iter().enumerate() { ui.label(i.to_string()); for value in row { ui.label(value); } ui.end_row(); }
                            });
                        });
                        ui.small(if table.rsa.is_some() { "RSA rows are shown in ROI-ID order; the source file is unchanged." } else { "Cells in this inspector are limited to 160 characters; the source file is unchanged." });
                    });
                }
                if self.busy { ui.horizontal(|ui| { ui.spinner(); ui.label("Working in the background…"); }); }
                if let Some(error) = &self.error { ui.colored_label(Color32::from_rgb(170,55,45), error); }
                });
                ui.separator();
                ui.horizontal(|ui| {
                    let response = ui.add_enabled(!self.busy && self.preview.is_some() && !changed, egui::Button::new(RichText::new("Add graph to figure").color(Color32::WHITE)).fill(Color32::from_rgb(20,125,146)));
                    #[cfg(test)] { self.add_rect = Some(response.rect); }
                    add = response.clicked();
                    retry = ui.button("Reload file / Retry").clicked();
                    if !self.queue.is_empty() { skip = ui.button("Skip file").clicked(); }
                    close = ui.button("Cancel import").clicked();
                });
            });
        #[cfg(test)]
        {
            self.page_rect = Some(page.response.rect);
        }
        #[cfg(not(test))]
        let _ = page;
        if close {
            self.cancel();
            ctx.request_repaint();
        } else if skip {
            if let Some(worker) = &mut self.worker {
                worker.cancel();
            }
            self.advance();
        } else if retry {
            self.retry();
        } else if changed {
            self.refresh();
        } else if add && let Some(preview) = self.preview.take() {
            self.busy = true;
            self.worker.as_mut().unwrap().start(Work::Confirm(preview));
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn settle(import: &mut ImportUi, ctx: &egui::Context, size: [f32; 2]) {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        loop {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
                ..Default::default()
            };
            ctx.run_ui(input, |ui| {
                assert!(import.show(ui).is_none());
            })
            .drop_without_applying_deltas();
            if !import.busy {
                return;
            }
            assert!(std::time::Instant::now() < deadline, "preview timed out");
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
    }
    #[test]
    fn rsa_drop_review_selects_effects_and_filters_roi_ids() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/rsa/ordinary.rsa.1D");
        let ctx = egui::Context::default();
        let mut import = ImportUi::default();
        import.enqueue(vec![path], &ctx);
        settle(&mut import, &ctx, [1320., 900.]);
        let table = import.table.as_ref().unwrap();
        assert!(table.rsa.is_some());
        assert_eq!(table.numeric().len(), 2);
        assert_eq!(import.mapping.as_ref().unwrap().kind, Kind::Histogram);
        assert_eq!(import.preview.as_ref().unwrap().figure.plots.len(), 1);
        let mut mapping = Mapping::for_kind(table, Kind::Heat);
        mapping.roi_filter = "12,3".into();
        import.mapping = Some(mapping);
        import.refresh();
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.error.is_none(), "{:?}", import.error);
        let crate::model::Content::Heat {
            matrix, row_labels, ..
        } = &import.preview.as_ref().unwrap().figure.plots[0].content
        else {
            panic!()
        };
        assert_eq!(row_labels, &["12: Occipital", "3: Left temporal"]);
        assert_eq!(matrix[0], [0.18, 0.73]);
        import.mapping.as_mut().unwrap().roi_filter = "999".into();
        import.refresh();
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.preview.is_none());
        assert!(import.error.as_ref().unwrap().contains("absent"));
    }
    #[test]
    fn add_action_is_visible_and_clickable_without_scrolling_on_short_windows() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("data.1D");
        std::fs::write(&path, "1 2\n2 3\n3 4\n").unwrap();
        let ctx = egui::Context::default();
        let mut import = ImportUi::default();
        import.enqueue(vec![path], &ctx);
        settle(&mut import, &ctx, [850., 500.]);
        settle(&mut import, &ctx, [850., 500.]);
        assert_eq!(import.page_rect.unwrap(), ctx.content_rect());
        settle(&mut import, &ctx, [1440., 900.]);
        assert_eq!(import.page_rect.unwrap(), ctx.content_rect());
        settle(&mut import, &ctx, [850., 500.]);
        let action = import.add_rect.unwrap();
        assert!(
            ctx.content_rect().contains_rect(action),
            "Add button outside viewport: {action:?}"
        );
        let position = action.center();
        let mut accepted = None;
        for pressed in [true, false] {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    [850., 500.].into(),
                )),
                events: vec![
                    egui::Event::PointerMoved(position),
                    egui::Event::PointerButton {
                        pos: position,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                ..Default::default()
            };
            ctx.run_ui(input, |ui| {
                accepted = import.show(ui);
            })
            .drop_without_applying_deltas();
        }
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        while accepted.is_none() {
            assert!(
                std::time::Instant::now() < deadline,
                "click did not confirm preview"
            );
            std::thread::sleep(std::time::Duration::from_millis(5));
            ctx.run_ui(egui::RawInput::default(), |ui| {
                accepted = import.show(ui);
            })
            .drop_without_applying_deltas();
        }
        assert_eq!(accepted.unwrap().figure.plots.len(), 1);
    }
    #[test]
    fn queued_imports_skip_retry_validate_changes_and_cancel_without_adding() {
        let d = tempfile::tempdir().unwrap();
        let bad = d.path().join("bad.1D");
        let good = d.path().join("good.tsv");
        std::fs::write(&bad, "1 2\n3\n").unwrap();
        std::fs::write(&good, "label\tx\ty\nA\t0\t1\nB\t1\t3\n").unwrap();
        let ctx = egui::Context::default();
        let mut import = ImportUi::default();
        import.enqueue(vec![bad.clone(), good.clone()], &ctx);
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.error.is_some());
        assert_eq!(import.queue.len(), 1);
        import.worker.as_mut().unwrap().cancel();
        import.advance();
        settle(&mut import, &ctx, [1320., 900.]);
        assert!(import.preview.is_some());
        assert_eq!(import.table.as_ref().unwrap().numeric(), [1, 2]);
        // The file can change after preview. Confirmation must refuse stale data.
        std::fs::write(&good, "label\tx\ty\nA\t0\t2\nB\t1\t4\n").unwrap();
        let prepared = import.preview.take().unwrap();
        import.busy = true;
        import
            .worker
            .as_mut()
            .unwrap()
            .start(Work::Confirm(prepared));
        settle(&mut import, &ctx, [1320., 900.]);
        assert!(import.error.as_ref().unwrap().contains("changed"));
        import.retry();
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.error.is_none());
        assert!(import.preview.is_some());
        import.cancel();
        assert!(!import.visible());
        assert!(import.preview.is_none());
        assert!(import.queue.is_empty());
        import.enqueue(vec![bad; 17], &ctx);
        assert!(!import.visible());
        assert!(import.error.as_ref().unwrap().contains("16"));
    }
}
