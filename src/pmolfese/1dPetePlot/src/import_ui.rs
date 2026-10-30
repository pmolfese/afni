//! Native-window import workflow. Temporary candidates never mutate the document.
use crate::import::{Answer, Kind, Mapping, Prepared, Table, Work, Worker};
use eframe::egui::{self, Color32, RichText};
use std::{collections::VecDeque, path::PathBuf, sync::Arc};

#[derive(Default)]
pub struct ImportUi {
    pub edit_target: Option<usize>,
    edit_recipe: Option<(
        crate::cli::Options,
        crate::document::Roles,
        crate::document::Panel,
    )>,
    worker: Option<Worker>,
    queue: VecDeque<PathBuf>,
    path: Option<PathBuf>,
    table: Option<Arc<Table>>,
    mapping: Option<Mapping>,
    retry_mapping: Option<(PathBuf, Vec<String>, Mapping)>,
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
    pub fn edit(
        &mut self,
        id: usize,
        o: crate::cli::Options,
        roles: crate::document::Roles,
        panel: crate::document::Panel,
        ctx: &egui::Context,
    ) -> anyhow::Result<()> {
        let file = o
            .files
            .first()
            .ok_or_else(|| anyhow::anyhow!("choose an ordinary table graph"))?;
        let path = PathBuf::from(crate::document::split_spec(file).0);
        self.cancel();
        self.edit_target = Some(id);
        self.edit_recipe = Some((o, roles, panel));
        self.enqueue(vec![path], ctx);
        Ok(())
    }
    pub fn cancel(&mut self) {
        if let Some(w) = &mut self.worker {
            w.cancel();
        }
        self.edit_target = None;
        self.edit_recipe = None;
        self.retry_mapping = None;
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
        self.retry_mapping =
            self.table
                .as_ref()
                .zip(self.mapping.as_ref())
                .map(|(table, mapping)| {
                    (
                        table.path.clone(),
                        table.columns.iter().map(|c| c.name.clone()).collect(),
                        mapping.clone(),
                    )
                });
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
            worker.start(Work::Preview(table.clone(), Box::new(mapping.clone())));
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
                    let names: Vec<_> = table.columns.iter().map(|c| c.name.clone()).collect();
                    self.mapping = Some(
                        if let Some((path, old_names, mapping)) = self.retry_mapping.take()
                            && path == table.path
                            && old_names == names
                        {
                            mapping
                        } else if let Some((o, roles, panel)) = &self.edit_recipe {
                            match Mapping::restore(&table, o, roles, panel.clone()) {
                                Ok(m) => m,
                                Err(e) => {
                                    self.error = Some(format!("{e:#}"));
                                    self.table = Some(table);
                                    return None;
                                }
                            }
                        } else {
                            Mapping::for_kind(
                                &table,
                                if let Some(long) = table.long() {
                                    if long.time { Kind::Line } else { Kind::Scatter }
                                } else if table.rsa.is_some() {
                                    Kind::Histogram
                                } else {
                                    Kind::Line
                                },
                            )
                        },
                    );
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
                ui.label(RichText::new(if self.edit_target.is_some(){"EDIT SCIENTIFIC ATTACHMENTS"}else{"IMPORT DATA"}).size(11.).color(Color32::from_rgb(20,125,146)));
                ui.heading(&filename);
                ui.label(if self.edit_target.is_some(){"Edit the scientific attachments, preview, then Apply. Existing titles and axis labels are retained; review their meaning. Cancel keeps the current graph."}else{"Preview your graph, then add it to the figure. Cancel returns to your figure."});
                ui.small(path.display().to_string());
                if !self.queue.is_empty() { ui.label(format!("{} more file(s) queued · each gets its own preview", self.queue.len())); }
                if let (Some(table), Some(mapping)) = (&self.table, &mut self.mapping) {
                    ui.label(format!("{} rows · {} columns · {} {}", table.rows, table.columns.len(), table.numeric().len(), if table.rsa.is_some() { "RSA effects" } else { "numeric" }));
                    if table.rsa.is_some() {
                        ui.label(if table.long().is_some() { "Long 3dRSA results · identities retained; statistics are metadata, never effects." } else { "3dRSA results · choose effects of the same kind (r, b, or pr). Statistics are retained as metadata." });
                        ui.horizontal(|ui| { ui.label(if table.long().is_some() { "ROI ID (one per graph)" } else { "ROI IDs (blank = all)" }); changed |= ui.text_edit_singleline(&mut mapping.roi_filter).on_hover_text("Comma-separated ROI IDs, e.g. 3,12,27. The chosen order is retained.").changed(); });
                    }
                    if let Some(rsa)=&table.rsa && rsa.long.is_none() {
                        let previous=mapping.rsa_family.clone();
                        egui::ComboBox::from_id_salt("rsa-family").selected_text(crate::rsa_extra::FAMILIES.iter().find(|(f,_)|*f==mapping.rsa_family).map(|(_,l)|*l).unwrap_or("Primary effects")).show_ui(ui,|ui|{
                            for (family,label) in crate::rsa_extra::FAMILIES {if family=="effects" || (0..rsa.columns.len()).any(|i|crate::rsa_extra::family(rsa,i)==Some(family)){ui.selectable_value(&mut mapping.rsa_family,family.into(),label);}}
                        });
                        if previous!=mapping.rsa_family {let family=mapping.rsa_family.clone();mapping.set_rsa_family(table,&family);changed=true;}
                    }
                    if let Some(long) = table.long() {
                        if !long.time { ui.horizontal(|ui| { ui.label("Models (blank = all)"); changed |= ui.text_edit_singleline(&mut mapping.rsa_models).changed(); });
                            ui.small(format!("Available models: {}", table.rsa.as_ref().unwrap().models.join(", "))); }
                        ui.horizontal(|ui| { ui.label(if long.time { "Time indices (blank = all)" } else { "Summaries (blank = all)" }); changed |= ui.text_edit_singleline(&mut mapping.rsa_steps).changed(); });
                        ui.small(format!("Available: {}{}", long.steps.iter().take(24).cloned().collect::<Vec<_>>().join(", "), if long.steps.len()>24 { " …" } else { "" }));
                    }
                    if self.edit_target.is_some() {changed |= ui.checkbox(&mut mapping.reset_style,"Reset styling for changed series roles").changed();}
                    ui.horizontal_wrapped(|ui| {
                        for kind in Kind::ALL {
                            let enabled = if mapping.rsa_family!="effects" {match kind {Kind::Histogram|Kind::Heat=>true,Kind::Intervals=>["commonality","partial-r2"].contains(&mapping.rsa_family.as_str()),Kind::Scatter=>table.rsa.as_ref().is_some_and(|t|(0..t.columns.len()).filter(|i|crate::rsa_extra::family(t,*i)==Some(mapping.rsa_family.as_str())).count()>=2),_=>false}} else if let Some(long) = table.long() { match kind { Kind::Line | Kind::Intervals => long.time, Kind::Scatter => !long.time, Kind::Heat => true, _ => false } } else { match kind { Kind::Rdm => table.rsa.is_none() && table.square() && table.path.extension().is_some_and(|e| e.eq_ignore_ascii_case("1d")), Kind::Intervals => table.rsa.is_some(), Kind::Line if table.rsa.is_some() => false, Kind::Scatter => table.rsa.as_ref().map(|rsa| ["r", "b", "pr"].iter().any(|kind| rsa.effects.keys().filter(|(_, k)| k == kind).count() >= 2)).unwrap_or_else(|| table.numeric().len() >= 2), Kind::Heat => table.rsa.is_some() || table.square(), _ => !table.numeric().is_empty() } };
                            if ui.add_enabled(enabled, egui::Button::new(kind.name()).selected(mapping.kind == kind)).on_hover_text(kind.explanation()).clicked() {
                                let missing = mapping.missing.clone(); let rois = mapping.roi_filter.clone(); let ci = mapping.ci.clone(); let rdm = mapping.rdm.clone(); let attachments = mapping.attachments.clone(); let models = mapping.rsa_models.clone(); let steps = mapping.rsa_steps.clone(); let band = mapping.rsa_band; let presentation=mapping.presentation.clone(); let reset=mapping.reset_style; let family=mapping.rsa_family.clone(); *mapping = Mapping::for_kind(table, kind); mapping.missing = missing; mapping.roi_filter = rois; mapping.ci = ci; mapping.rdm = rdm; mapping.attachments = attachments; mapping.rsa_models = models; mapping.rsa_steps = steps; mapping.rsa_band = band; mapping.presentation=presentation; mapping.reset_style=reset; if family!="effects"{mapping.set_rsa_family(table,&family);} changed = true;
                            }
                        }
                    });
                    let explanation = if mapping.rsa_family!="effects" { "Named quantities in one scientific family. Negative commonality is retained; weights are descriptive. No refitting, significance filtering, normalization or inferred intervals." } else if table.long().is_some() { "One ROI per graph. Time indices retain source order; runs, MEAN and contrasts are categorical. Full identities and supplied statistics appear on hover." } else if table.rsa.is_some() { match mapping.kind {
                        Kind::Intervals => "One model with supplied CI endpoints, up to 40 ROIs. Beta intervals cannot be applied to partial r.",
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
                            if table.long().is_none() && matches!(mapping.kind, Kind::Line | Kind::Scatter) {
                                ui.label("X coordinates");
                                egui::ComboBox::from_id_salt("import-x").selected_text(mapping.x.map(|i| table.columns[i].name.as_str()).unwrap_or("Row index (0-based)"))
                                    .width(225.).show_ui(ui, |ui| {
                                        if mapping.kind == Kind::Line { changed |= ui.selectable_value(&mut mapping.x, None, "Row index (0-based)").changed(); }
                                        for i in table.numeric() { if !mapping.allows(table,i) {continue;} changed |= ui.selectable_value(&mut mapping.x, Some(i), &table.columns[i].name).changed(); }
                                    });
                            }
                            if table.long().is_some_and(|l| l.time) && mapping.kind == Kind::Line {
                                changed |= ui.checkbox(&mut mapping.rsa_band, "Show supplied CI band").changed();
                            }
                            if mapping.kind == Kind::Intervals || table.long().is_some() && mapping.kind == Kind::Line && mapping.rsa_band {
                                ui.label("Confidence interval method");
                                egui::ComboBox::from_id_salt("import-ci").selected_text(match mapping.ci.as_str() { "boot" => "Subject bootstrap", "cboot" => "Condition bootstrap", "dual" => "Dual bootstrap", _ => "Auto (one available method)" })
                                    .width(225.).show_ui(ui, |ui| {
                                        for (value, label) in [("auto", "Auto (one available method)"), ("boot", "Subject bootstrap"), ("cboot", "Condition bootstrap"), ("dual", "Dual bootstrap")] {
                                            changed |= ui.selectable_value(&mut mapping.ci, value.into(), label).changed();
                                        }
                                    });
                                ui.small("Gap policy retains incomplete rows marked missing, without a point or interval.");
                            }
                            if mapping.kind == Kind::Rdm {
                                for (label, value) in [("Condition labels (required)", &mut mapping.rdm.labels), ("Display order (optional)", &mut mapping.rdm.order)] {
                                    ui.label(label);
                                    changed |= ui.text_edit_singleline(value).changed();
                                    if ui.button(format!("Choose {label}…")).clicked()
                                        && let Some(path) = rfd::FileDialog::new().set_title(label).pick_file() {
                                        *value = path.to_string_lossy().into(); changed = true;
                                    }
                                }
                                ui.small("One unique label per line, in source matrix order. Display order must contain the same names.");
                                egui::ComboBox::from_id_salt("rdm-sense").selected_text(&mapping.rdm.sense).show_ui(ui, |ui| {
                                    for sense in ["unspecified", "similarity", "dissimilarity"] { changed |= ui.selectable_value(&mut mapping.rdm.sense, sense.into(), sense).changed(); }
                                });
                                ui.label("Fixed color limits (optional MIN:MAX)");
                                changed |= ui.text_edit_singleline(&mut mapping.rdm.clim).changed();
                                ui.small("Use identical fixed limits for comparable graphs. Blank scales this matrix independently; automatic shared scales across files are available with CLI -rdm-scale shared.");
                                ui.label("The entire matrix is used; both axes follow condition identity.");
                            } else {
                            ui.label(if table.rsa.is_some() { "Y / model effect columns" } else if mapping.kind == Kind::Heat { "Matrix columns (row order retained)" } else { "Y / value columns" });
                            egui::ScrollArea::vertical().id_salt("import-columns").max_height(190.).show(ui, |ui| {
                                for (i, column) in table.columns.iter().enumerate() {
                                    let mut selected = mapping.y.contains(&i);
                                    let label = if !column.role.is_empty() { format!("{} · {}", column.name, column.role) } else { format!("{}{}", column.name, if column.invalid { " · infinity" } else if !column.numeric { " · text" } else { "" }) };
                                    let response = ui.add_enabled(column.numeric && !column.invalid && mapping.allows(table,i), egui::Checkbox::new(&mut selected, label));
                                    if response.changed() {
                                        if selected { mapping.y.push(i); mapping.y.sort_unstable(); } else { mapping.y.retain(|&c| c != i); } changed = true;
                                    }
                                    if column.missing > 0 { ui.small(format!("{} missing", column.missing)); }
                                }
                            });
                            }
                            changed |= mapping.attachments.controls(ui, table, mapping.kind, &mapping.y);
                            ui.separator(); ui.label("Missing observations");
                            egui::ComboBox::from_id_salt("import-missing").selected_text(match mapping.missing.as_str() { "gap" => "Allow gaps / omit from counts", "mean" => "Impute column mean", _ => "Reject missing values" })
                                .width(225.).show_ui(ui, |ui| {
                                    for (value,label) in [("error","Reject missing values"),("gap","Allow gaps / omit from counts"),("mean","Impute column mean")] {
                                        if (table.rsa.is_some() || mapping.kind == Kind::Rdm || mapping.attachments.group.is_some() || matches!(mapping.attachments.uncertainty.as_str(), "widths" | "bounds" | "file")) && value == "mean" { continue; }
                                        changed |= ui.selectable_value(&mut mapping.missing, value.into(), label).changed();
                                    }
                                });
                            ui.small(if mapping.kind == Kind::Rdm { "Missing mirrored cells must match. No imputation, diagonal replacement, or similarity conversion." } else if table.rsa.is_some() { "Missing scatter pairs are omitted together with gap policy. No significance filtering or imputation." } else { "Text columns stay in the table; they are not converted to numbers. X must be finite." });
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
                        ui.small(if table.long().is_some() { "Long-table inspector retains source row order; plots use explicit identities." } else if table.rsa.is_some() { "RSA rows are shown in ROI-ID order; the source file is unchanged." } else { "Cells in this inspector are limited to 160 characters; the source file is unchanged." });
                    });
                }
                if self.busy { ui.horizontal(|ui| { ui.spinner(); ui.label("Working in the background…"); }); }
                if let Some(error) = &self.error { ui.colored_label(Color32::from_rgb(170,55,45), error); }
                });
                ui.separator();
                ui.horizontal(|ui| {
                    let response = ui.add_enabled(!self.busy && self.preview.is_some() && !changed, egui::Button::new(RichText::new(if self.edit_target.is_some(){"Apply changes to graph"}else{"Add graph to figure"}).color(Color32::WHITE)).fill(Color32::from_rgb(20,125,146)));
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
    fn editing_attachment_review_restores_roles_and_requires_alignment() {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/attachments");
        let table = crate::import::inspect(&base.join("signals.tsv"), || false).unwrap();
        let mut m = Mapping::for_kind(&table, Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.censor_file = base.join("censor.1D").to_string_lossy().into();
        m.attachments.aligned = true;
        let p = crate::import::prepare(&table, &m, || false).unwrap();
        let o = p.origin.verify().unwrap();
        let panel = crate::document::Presentation::capture(&p.figure).panels[0].clone();
        let ctx = egui::Context::default();
        let mut editor = ImportUi::default();
        editor.edit(7, o, p.roles, panel, &ctx).unwrap();
        settle(&mut editor, &ctx, [850., 500.]);
        assert_eq!(editor.edit_target, Some(7));
        assert!(editor.preview.is_none());
        assert!(editor.error.as_ref().unwrap().contains("confirm"));
        let mapping = editor.mapping.as_mut().unwrap();
        assert_eq!(mapping.x, Some(0));
        assert_eq!(mapping.y, [1]);
        assert!(mapping.attachments.censor_file.ends_with("censor.1D"));
        mapping.attachments.aligned = true;
        mapping.attachments.uncertainty = "widths".into();
        mapping.attachments.bands.insert(1, [Some(4), None]);
        editor.refresh();
        settle(&mut editor, &ctx, [850., 500.]);
        assert!(editor.preview.is_some(), "{:?}", editor.error);
        assert!(ctx.content_rect().contains_rect(editor.add_rect.unwrap()));
        editor.cancel();
        assert!(editor.edit_target.is_none());
    }
    #[test]
    fn fit_and_commonality_review_selects_only_the_chosen_scientific_family() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("examples/rsa/fitted-commonality.rsa.1D");
        let ctx = egui::Context::default();
        let mut editor = ImportUi::default();
        editor.enqueue(vec![path], &ctx);
        settle(&mut editor, &ctx, [850., 600.]);
        for (family, n) in [
            ("fit", 2),
            ("weights", 2),
            ("commonality", 3),
            ("partial-r2", 2),
        ] {
            let table = editor.table.as_ref().unwrap();
            let mut m = Mapping::for_kind(table, Kind::Heat);
            m.set_rsa_family(table, family);
            assert_eq!(m.y.len(), n);
            editor.mapping = Some(m);
            editor.refresh();
            settle(&mut editor, &ctx, [850., 600.]);
            assert!(editor.preview.is_some(), "{:?}", editor.error);
        }
        let table = editor.table.as_ref().unwrap();
        let mut m = Mapping::for_kind(table, Kind::Intervals);
        m.set_rsa_family(table, "commonality");
        m.y = vec![
            table
                .columns
                .iter()
                .position(|c| c.name == "common_A_B")
                .unwrap(),
        ];
        editor.mapping = Some(m);
        editor.refresh();
        settle(&mut editor, &ctx, [850., 500.]);
        assert!(editor.preview.is_some());
        let crate::model::Content::Intervals { rows, .. } =
            &editor.preview.as_ref().unwrap().figure.plots[0].content
        else {
            panic!()
        };
        assert!(rows.iter().all(|r| r.values.unwrap()[0] < 0.));
    }
    #[test]
    fn long_rsa_review_defaults_and_filters_keep_roles_and_visible_add() {
        for (name, kind) in [("series", Kind::Line), ("runs-joint", Kind::Scatter)] {
            let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join(format!("examples/rsa/{name}.rsa.1D"));
            let ctx = egui::Context::default();
            let mut import = ImportUi::default();
            import.enqueue(vec![path], &ctx);
            settle(&mut import, &ctx, [850., 500.]);
            assert!(import.error.is_none(), "{:?}", import.error);
            assert_eq!(import.mapping.as_ref().unwrap().kind, kind);
            assert_eq!(import.preview.as_ref().unwrap().figure.plots.len(), 1);
            assert!(ctx.content_rect().contains_rect(import.add_rect.unwrap()));
            let mapping = import.mapping.as_mut().unwrap();
            mapping.roi_filter = "12".into();
            if kind == Kind::Line {
                mapping.rsa_steps = "10,2,0".into();
                mapping.rsa_band = true;
            } else {
                mapping.rsa_models = "visual".into();
                mapping.rsa_steps = "task,rest".into();
            }
            import.refresh();
            settle(&mut import, &ctx, [850., 500.]);
            assert!(import.error.is_none(), "{:?}", import.error);
            let plot = &import.preview.as_ref().unwrap().figure.plots[0];
            assert!(plot.title.starts_with("12:"));
            if let crate::model::Content::Lines(ss) = &plot.content {
                assert_eq!(ss[0].points.len(), 3);
                assert!(ss[0].band.is_some());
            } else if let crate::model::Content::Scatter(ss) = &plot.content {
                assert_eq!(ss.len(), 1);
                assert_eq!(ss[0].label, "visual");
                assert_eq!(ss[0].points.len(), 2);
            } else {
                panic!()
            }
            import.mapping.as_mut().unwrap().roi_filter = "3,12".into();
            import.refresh();
            settle(&mut import, &ctx, [850., 500.]);
            assert!(import.preview.is_none());
            assert!(import.error.as_ref().unwrap().contains("one ROI"));
        }
    }
    #[test]
    fn rdm_review_requires_labels_and_previews_order_and_fixed_scale() {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/rdm");
        let ctx = egui::Context::default();
        let mut import = ImportUi::default();
        import.enqueue(vec![base.join("model.1D")], &ctx);
        settle(&mut import, &ctx, [850., 600.]);
        import.mapping = Some(Mapping::for_kind(import.table.as_ref().unwrap(), Kind::Rdm));
        import.refresh();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.preview.is_none());
        assert!(
            import
                .error
                .as_ref()
                .unwrap()
                .contains("condition-label file")
        );
        let m = import.mapping.as_mut().unwrap();
        m.rdm.labels = base.join("conditions.txt").to_str().unwrap().into();
        m.rdm.order = base.join("display-order.txt").to_str().unwrap().into();
        m.rdm.clim = "-1:1".into();
        import.refresh();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.error.is_none(), "{:?}", import.error);
        let crate::model::Content::Heat {
            matrix,
            labels,
            clim,
            ..
        } = &import.preview.as_ref().unwrap().figure.plots[0].content
        else {
            panic!()
        };
        assert_eq!(matrix[0][1], 0.18);
        assert_eq!(&labels[..2], ["Faces", "Animals"]);
        assert_eq!(*clim, [-1., 1.]);
        assert!(import.add_rect.is_some_and(|r| r.bottom() <= 500.));
        import.mapping.as_mut().unwrap().rdm.clim = "1:-1".into();
        import.refresh();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.preview.is_none());
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
        assert_eq!(table.rsa.as_ref().unwrap().effects.len(), 2);
        assert!(table.numeric().len() > 2); // Fit/weight families are now available explicitly.
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
        let mut mapping = Mapping::for_kind(import.table.as_ref().unwrap(), Kind::Intervals);
        mapping.roi_filter = "12,3".into();
        import.mapping = Some(mapping);
        import.refresh();
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.error.is_none(), "{:?}", import.error);
        let crate::model::Content::Intervals { rows, method, .. } =
            &import.preview.as_ref().unwrap().figure.plots[0].content
        else {
            panic!()
        };
        assert_eq!(rows[0].label, "12: Occipital");
        assert_eq!(rows[0].values, Some([0.18, -0.04, 0.34]));
        assert!(method.contains("95%"));
        import.mapping.as_mut().unwrap().ci = "dual".into();
        import.refresh();
        settle(&mut import, &ctx, [850., 600.]);
        assert!(import.preview.is_none());
        assert!(import.error.as_ref().unwrap().contains("no bounds"));
        import.mapping.as_mut().unwrap().ci = "auto".into();
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
    #[test]
    fn scientific_attachment_preview_requires_alignment_and_rejects_stale_confirmation() {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("examples/attachments");
        let d = tempfile::tempdir().unwrap();
        let censor = d.path().join("mask.1D");
        std::fs::copy(base.join("censor.1D"), &censor).unwrap();
        let ctx = egui::Context::default();
        let mut import = ImportUi::default();
        import.enqueue(vec![base.join("signals.tsv")], &ctx);
        settle(&mut import, &ctx, [850., 500.]);
        let m = import.mapping.as_mut().unwrap();
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.uncertainty = "bounds".into();
        m.attachments.bands.insert(1, [Some(2), Some(3)]);
        m.attachments.censor_file = censor.display().to_string();
        m.attachments.xmat = base.join("X.mat.1D").display().to_string();
        import.refresh();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.preview.is_none());
        assert!(import.error.as_ref().unwrap().contains("confirm"));
        import.mapping.as_mut().unwrap().attachments.aligned = true;
        import.refresh();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.error.is_none(), "{:?}", import.error);
        assert!(import.preview.is_some());
        assert!(ctx.content_rect().contains_rect(import.add_rect.unwrap()));
        std::fs::write(&censor, "1\n1\n1\n1\n1\n1\n").unwrap();
        let position = import.add_rect.unwrap().center();
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
                        modifiers: Default::default(),
                    },
                ],
                ..Default::default()
            };
            ctx.run_ui(input, |ui| assert!(import.show(ui).is_none()))
                .drop_without_applying_deltas();
        }
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.error.as_ref().unwrap().contains("source changed"));
        assert!(import.preview.is_none());
        import.retry();
        settle(&mut import, &ctx, [850., 500.]);
        assert!(import.preview.is_some(), "{:?}", import.error);
        let m = import.mapping.as_ref().unwrap();
        assert_eq!(m.y, vec![1]);
        assert_eq!(m.attachments.uncertainty, "bounds");
        assert_eq!(m.attachments.censor_file, censor.display().to_string());
        import.cancel();
        assert!(!import.visible());
    }
}
