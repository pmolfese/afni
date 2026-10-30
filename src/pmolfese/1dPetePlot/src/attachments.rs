//! Scientific roles for the graph being imported. Paths live in its fingerprinted recipe.
use crate::import::{Kind, Table};
use anyhow::{Context, Result, ensure};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Attachments {
    pub group: Option<usize>,
    /// none, widths, bounds, or file; supplied quantities, never estimated here.
    pub uncertainty: String,
    /// Original Y column -> width or (lower, upper) original columns.
    pub bands: BTreeMap<usize, [Option<usize>; 2]>,
    pub width_file: String,
    pub include_bounds: bool,
    pub censor_file: String,
    pub censor_column: Option<usize>,
    pub xmat: String,
    pub threshold_percent: u32,
    pub palette: String,
    pub aligned: bool,
}
impl Attachments {
    pub fn has_sidecars(&self) -> bool {
        (!self.width_file.trim().is_empty() && self.uncertainty == "file")
            || !self.censor_file.trim().is_empty()
            || !self.xmat.trim().is_empty()
    }
    pub fn active(&self) -> bool {
        self.group.is_some()
            || (!self.uncertainty.is_empty() && self.uncertainty != "none")
            || self.censor_column.is_some()
            || self.has_sidecars()
    }
    pub fn arguments(
        &self,
        table: &Table,
        kind: Kind,
        x: Option<usize>,
        y: &[usize],
        missing: &str,
    ) -> Result<Vec<String>> {
        if !self.active() {
            return Ok(vec![]);
        }
        ensure!(
            table.rsa.is_none() && matches!(kind, Kind::Line | Kind::Scatter | Kind::Histogram),
            "scientific attachments require an ordinary line, scatter, or distribution graph; clear attachments before changing plot family"
        );
        let bands = !self.uncertainty.is_empty() && self.uncertainty != "none";
        ensure!(
            !bands || kind == Kind::Line,
            "uncertainty bands require a line graph"
        );
        ensure!(
            kind == Kind::Line
                || (self.censor_column.is_none()
                    && self.censor_file.trim().is_empty()
                    && self.xmat.trim().is_empty()),
            "censor and X-matrix attachments require a line timeline"
        );
        ensure!(
            !self.has_sidecars() || self.aligned,
            "confirm that attachment rows use this graph's full sample order/timeline; matching lengths alone cannot establish alignment"
        );
        ensure!(
            !(bands || self.group.is_some()) || missing != "mean",
            "uncertainty/group roles do not support mean imputation; choose reject or gaps"
        );
        let file = table.path.to_str().context("path must be UTF-8")?;
        let mut args = vec![];
        let mut used = y.to_vec();
        used.extend(x);
        if let Some(g) = self.group {
            ensure!(
                g < table.columns.len() && !used.contains(&g),
                "group must be a distinct original table column"
            );
            used.push(g);
            args.extend(["-group-by".into(), g.to_string()]);
        }
        if let Some(c) = self.censor_column {
            ensure!(
                self.censor_file.trim().is_empty(),
                "choose a censor column or a censor file, not both"
            );
            ensure!(
                table
                    .columns
                    .get(c)
                    .is_some_and(|c| c.numeric && !c.invalid)
                    && !used.contains(&c),
                "censor must be a distinct numeric column"
            );
            used.push(c);
            args.extend(["-censor".into(), format!("{file}[{c}]")]);
        }
        if bands {
            match self.uncertainty.as_str() {
                "widths" | "bounds" => {
                    let mut first = vec![];
                    let mut second = vec![];
                    for &value in y {
                        let pair = self
                            .bands
                            .get(&value)
                            .context("assign uncertainty columns for every selected Y series")?;
                        for (slot, column) in pair
                            .iter()
                            .enumerate()
                            .take(if self.uncertainty == "bounds" { 2 } else { 1 })
                        {
                            let c = column.context("choose each supplied uncertainty column")?;
                            ensure!(
                                table
                                    .columns
                                    .get(c)
                                    .is_some_and(|c| c.numeric && !c.invalid)
                                    && !used.contains(&c),
                                "uncertainty columns must be distinct from X, Y, group, and censor"
                            );
                            if slot == 0 {
                                first.push(c.to_string());
                            } else {
                                second.push(c.to_string());
                            }
                        }
                    }
                    args.extend([
                        if self.uncertainty == "bounds" {
                            "-yfiles_lo"
                        } else {
                            "-yfiles_pm"
                        }
                        .into(),
                        format!("{file}[{}]", first.join(",")),
                    ]);
                    if self.uncertainty == "bounds" {
                        args.extend(["-yfiles_hi".into(), format!("{file}[{}]", second.join(","))]);
                    }
                }
                "file" => {
                    ensure!(
                        !self.width_file.trim().is_empty(),
                        "choose a supplied width file"
                    );
                    args.extend(["-yfiles_pm".into(), self.width_file.clone()]);
                }
                _ => anyhow::bail!("unknown uncertainty role"),
            }
            if self.include_bounds {
                args.push("-ylim_use_pm".into());
            }
        }
        if !self.censor_file.trim().is_empty() {
            args.extend(["-censor".into(), self.censor_file.clone()]);
        }
        if !self.xmat.trim().is_empty() {
            ensure!(
                self.threshold_percent < 100,
                "X-matrix threshold must be less than 100%"
            );
            args.extend([
                "-shade-xmat".into(),
                self.xmat.clone(),
                "-shade-xmat-threshold".into(),
                (self.threshold_percent as f64 / 100.).to_string(),
            ]);
            if !self.palette.trim().is_empty() {
                args.push("-shade-xmat-colors".into());
                args.extend(self.palette.split_whitespace().map(str::to_owned));
            }
        }
        Ok(args)
    }
    /// Check bounded local sidecars before the existing scientific readers allocate data.
    pub fn validate_files(&self) -> Result<()> {
        for path in [
            if self.uncertainty == "file" {
                &self.width_file
            } else {
                ""
            },
            &self.censor_file,
            &self.xmat,
        ] {
            if path.trim().is_empty() {
                continue;
            }
            let p = std::path::Path::new(path);
            ensure!(
                !path.ends_with([']', '}', '\'']) && !path.starts_with("1D:") && path != "-",
                "interactive sidecars need regular, unselected files; use table-column roles or prepare a separate file"
            );
            let meta =
                std::fs::metadata(p).with_context(|| format!("opening attachment {path}"))?;
            ensure!(
                meta.is_file() && meta.len() <= crate::import::MAX_BYTES,
                "attachment must be a regular file up to 16 MiB"
            );
        }
        Ok(())
    }
}

impl Attachments {
    pub fn controls(
        &mut self,
        ui: &mut eframe::egui::Ui,
        table: &Table,
        kind: Kind,
        y: &[usize],
    ) -> bool {
        use eframe::egui;
        let before = self.clone();
        ui.separator();
        ui.strong("Scientific roles & attachments");
        ui.small("Target: this preview → a new graph. Existing graphs are unchanged.");
        if self.active() && ui.button("Clear scientific roles / attachments").clicked() {
            *self = Self {
                include_bounds: true,
                threshold_percent: 5,
                ..Default::default()
            };
        }
        if table.rsa.is_some() || !matches!(kind, Kind::Line | Kind::Scatter | Kind::Histogram) {
            if self.active() {
                ui.label("Clear these roles before using this plot family.");
            }
            return *self != before;
        }
        column(ui, "Group column", &mut self.group, table, false);
        if self.group.is_some() {
            ui.small("Separate series by category, in first-seen order. No averaging or interpolation. Lines break between groups; missing groups are rejected. Up to 16 total series.");
        }
        if kind == Kind::Line {
            ui.label("Supplied uncertainty");
            egui::ComboBox::from_id_salt("uncertainty-mode")
                .selected_text(match self.uncertainty.as_str() {
                    "widths" => "± width columns",
                    "bounds" => "Lower / upper columns",
                    "file" => "± width file",
                    _ => "None",
                })
                .show_ui(ui, |ui| {
                    for (value, label) in [
                        ("none", "None"),
                        ("widths", "± width columns"),
                        ("bounds", "Lower / upper columns"),
                        ("file", "± width file"),
                    ] {
                        ui.selectable_value(&mut self.uncertainty, value.into(), label);
                    }
                });
            if matches!(self.uncertainty.as_str(), "widths" | "bounds") {
                for &value in y {
                    ui.push_id(value, |ui| {
                        ui.label(format!("For Y: {}", table.columns[value].name));
                        let role = self.bands.entry(value).or_insert([None, None]);
                        column(
                            ui,
                            if self.uncertainty == "bounds" {
                                "Lower endpoint"
                            } else {
                                "Absolute ± width"
                            },
                            &mut role[0],
                            table,
                            true,
                        );
                        if self.uncertainty == "bounds" {
                            column(ui, "Upper endpoint", &mut role[1], table, true);
                        }
                    });
                }
            } else if self.uncertainty == "file" {
                file(ui, "Uncertainty widths", &mut self.width_file);
                ui.small("One width column per selected Y, in Y selection order; full row count required. Widths use absolute values.");
            }
            if !self.uncertainty.is_empty() && self.uncertainty != "none" {
                ui.checkbox(
                    &mut self.include_bounds,
                    "Include bands in automatic Y range",
                );
                ui.small("Supplied values only: no SE/CI calculation or assumed confidence level. Missing bounds break the band; observed estimates remain visible. No mean imputation.");
            }
            ui.separator();
            column(ui, "Censor column", &mut self.censor_column, table, true);
            file(ui, "Censor mask", &mut self.censor_file);
            ui.small("Choose column OR file: exactly 0 = censored, 1 = retained. Shading preserves every observation. X must increase.");
            file(ui, "AFNI X matrix", &mut self.xmat);
            if !self.xmat.trim().is_empty() {
                ui.horizontal(|ui| {
                    ui.label("Activity threshold (%)");
                    ui.add(egui::DragValue::new(&mut self.threshold_percent).range(0..=99));
                });
                ui.label("Condition colors (optional, space separated)");
                ui.text_edit_singleline(&mut self.palette).on_hover_text("For example: #147d92 #e18a47 #7064b5. Positive condition groups follow numeric group order.");
                ui.small("Modeled response above the chosen fraction of each condition's peak magnitude—not stimulus onset timing. Baseline/motion groups are excluded. NRowFull / GoodList must match; missing matrix rows remain unshaded.");
            }
        } else if self.has_sidecars()
            || self.censor_column.is_some()
            || (!self.uncertainty.is_empty() && self.uncertainty != "none")
        {
            ui.label(
                "Bands and timeline attachments need Lines. Clear roles or switch back to Lines.",
            );
        }
        if self.width_file != before.width_file
            || self.censor_file != before.censor_file
            || self.xmat != before.xmat
        {
            self.aligned = false;
        }
        if self.has_sidecars() {
            ui.add_space(8.);
            ui.checkbox(
                &mut self.aligned,
                "I confirm these files use this graph's full sample order / timeline",
            );
            ui.small("No resampling, row dropping, or inferred joins. Matching lengths cannot prove that files belong together.");
        }
        *self != before
    }
}
fn column(
    ui: &mut eframe::egui::Ui,
    label: &str,
    selected: &mut Option<usize>,
    table: &Table,
    numeric: bool,
) {
    use eframe::egui;
    ui.label(label);
    egui::ComboBox::from_id_salt(label)
        .width(225.)
        .selected_text(
            selected
                .and_then(|i| table.columns.get(i))
                .map_or("None / choose…", |c| c.name.as_str()),
        )
        .show_ui(ui, |ui| {
            ui.selectable_value(selected, None, "None / choose…");
            for (i, c) in table.columns.iter().enumerate() {
                if !numeric || (c.numeric && !c.invalid) {
                    ui.selectable_value(selected, Some(i), &c.name);
                }
            }
        });
}
fn file(ui: &mut eframe::egui::Ui, label: &str, path: &mut String) {
    ui.label(label);
    ui.text_edit_singleline(path);
    ui.horizontal(|ui| {
        if ui.button(format!("Choose {label}…")).clicked()
            && let Some(p) = rfd::FileDialog::new()
                .set_title(label)
                .add_filter("AFNI tables", &["1D", "1d", "tsv"])
                .pick_file()
        {
            *path = p.to_string_lossy().into();
        }
        if !path.is_empty() && ui.small_button("Clear").clicked() {
            path.clear();
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        document,
        import::{self, Mapping},
        model::Content,
    };
    fn base() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/attachments")
    }
    fn mapping(table: &Table) -> Mapping {
        let mut m = Mapping::for_kind(table, Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.uncertainty = "bounds".into();
        m.attachments.bands.insert(1, [Some(2), Some(3)]);
        m.attachments.group = Some(5);
        m.attachments.censor_column = Some(6);
        m.attachments.xmat = base().join("X.mat.1D").display().to_string();
        m.attachments.aligned = true;
        m
    }
    #[test]
    fn combined_roles_keep_row_identity_asymmetric_bounds_censor_and_condition_colors() {
        let table = import::inspect(&base().join("signals.tsv"), || false).unwrap();
        let m = mapping(&table);
        let prepared = import::prepare(&table, &m, || false).unwrap();
        let plot = &prepared.figure.plots[0];
        let Content::Lines(series) = &plot.content else {
            panic!()
        };
        assert_eq!(series.len(), 2);
        assert_eq!(series[0].label, "effect · Rest");
        assert!(series[0].points[2][1].is_nan());
        assert_eq!(series[1].points[2], [4., 0.8]); // no row compression or censor deletion
        assert_eq!(series[1].band.as_ref().unwrap()[2], [0.6, 1.1]);
        assert!(series[0].band.as_ref().unwrap()[2][0].is_nan());
        assert_eq!(plot.qc.censored, [false, false, true, false, false, true]);
        assert_eq!(plot.qc.annotations.len(), 3); // gap at censored matrix row splits active conditions
        assert!(plot.qc.annotations.iter().all(|s| s.x[1] <= 9.));
        assert_ne!(
            plot.qc.annotations.first().unwrap().color,
            plot.qc.annotations.last().unwrap().color
        );
        assert!(plot.bounds.y[1] >= 1.1);
        let mut bad = m.clone();
        bad.attachments.aligned = false;
        assert!(
            import::prepare(&table, &bad, || false)
                .err()
                .unwrap()
                .to_string()
                .contains("confirm")
        );
        bad = m.clone();
        bad.missing = "mean".into();
        assert!(bad.arguments(&table).is_err());
        bad = m.clone();
        bad.attachments.bands.insert(1, [Some(0), Some(3)]);
        assert!(bad.arguments(&table).is_err());
        bad = m.clone();
        bad.kind = Kind::Scatter;
        assert!(bad.arguments(&table).is_err());
    }
    #[test]
    fn attachment_recipe_relocates_and_detects_changed_auxiliary_data_before_commit() {
        let dir = tempfile::tempdir().unwrap();
        let old = dir.path().join("old");
        std::fs::create_dir(&old).unwrap();
        for name in ["signals.tsv", "X.mat.1D", "censor.1D"] {
            std::fs::copy(base().join(name), old.join(name)).unwrap();
        }
        let table = import::inspect(&old.join("signals.tsv"), || false).unwrap();
        let mut m = mapping(&table);
        m.attachments.xmat = old.join("X.mat.1D").display().to_string();
        m.attachments.censor_column = None;
        m.attachments.censor_file = old.join("censor.1D").display().to_string();
        let prepared = import::prepare(&table, &m, || false).unwrap();
        let state = document::Snapshot {
            presentation: document::Presentation::capture(&prepared.figure),
            views: prepared.figure.plots.iter().map(|p| p.bounds).collect(),
            width: 900,
            height: 600,
            full_plot: false,
        };
        document::save(&old.join("plot.json"), &prepared.origin, &state, false).unwrap();
        let moved = dir.path().join("moved");
        std::fs::rename(old, &moved).unwrap();
        let loaded = document::load(&moved.join("plot.json")).unwrap();
        let renderer = crate::render::Renderer::new();
        assert_eq!(
            renderer
                .raster(
                    &renderer
                        .scene(&loaded.figure, &state.views, 900, 600)
                        .unwrap()
                )
                .unwrap()
                .data(),
            prepared.rgba
        );
        std::fs::write(moved.join("censor.1D"), "1\n1\n1\n1\n1\n1\n").unwrap();
        assert!(loaded.origin.verify().is_err()); // same check used by Confirm/Add
        assert!(document::load(&moved.join("plot.json")).is_err());
    }
    #[test]
    fn masks_bounds_and_design_reject_incompatible_or_missing_rows_without_imputing() {
        let dir = tempfile::tempdir().unwrap();
        let table = import::inspect(&base().join("signals.tsv"), || false).unwrap();
        for mask in ["1\n0\n", "1\n1\n2\n1\n1\n1\n", "1\n1\nNaN\n1\n1\n1\n"] {
            let file = dir.path().join("mask.1D");
            std::fs::write(&file, mask).unwrap();
            let mut m = mapping(&table);
            m.attachments.censor_column = None;
            m.attachments.censor_file = file.display().to_string();
            assert!(import::prepare(&table, &m, || false).is_err());
        }
        let mut m = mapping(&table);
        m.attachments.bands.insert(1, [Some(3), Some(2)]);
        assert!(import::prepare(&table, &m, || false).is_err());
        let file = dir.path().join("bounds.tsv");
        std::fs::write(
            &file,
            "x\ty\tlo\thi\n0\t1\t0\t2\n1\t2\tNaN\t3\n2\t3\t2\t5\n",
        )
        .unwrap();
        let table = import::inspect(&file, || false).unwrap();
        let mut m = Mapping::for_kind(&table, Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.missing = "gap".into();
        m.attachments.uncertainty = "bounds".into();
        m.attachments.bands.insert(1, [Some(2), Some(3)]);
        let prepared = import::prepare(&table, &m, || false).unwrap();
        let Content::Lines(s) = &prepared.figure.plots[0].content else {
            panic!()
        };
        assert_eq!(s[0].points[1], [1., 2.]);
        assert!(s[0].band.as_ref().unwrap()[1].iter().all(|x| x.is_nan()));
        m.attachments.xmat = base().join("X.mat.1D").display().to_string();
        m.attachments.aligned = true;
        assert!(
            import::prepare(&table, &m, || false)
                .err()
                .unwrap()
                .to_string()
                .contains("full 6-sample")
        );
    }
    #[test]
    fn width_columns_and_files_match_and_require_one_width_per_y() {
        let table = import::inspect(&base().join("signals.tsv"), || false).unwrap();
        let mut m = Mapping::for_kind(&table, Kind::Line);
        m.x = Some(0);
        m.y = vec![1];
        m.attachments.uncertainty = "widths".into();
        m.attachments.bands.insert(1, [Some(4), None]);
        let column = import::prepare(&table, &m, || false).unwrap();
        let d = tempfile::tempdir().unwrap();
        let widths = d.path().join("widths.1D");
        std::fs::write(&widths, "0.1\n-0.15\n0.2\n0.15\n0.2\n0.1\n").unwrap();
        m.attachments.uncertainty = "file".into();
        m.attachments.width_file = widths.display().to_string();
        m.attachments.aligned = true;
        let file = import::prepare(&table, &m, || false).unwrap();
        assert_eq!(file.rgba, column.rgba);
        std::fs::write(
            &widths,
            "0.1 0.2\n0.1 0.2\n0.1 0.2\n0.1 0.2\n0.1 0.2\n0.1 0.2\n",
        )
        .unwrap();
        assert!(import::prepare(&table, &m, || false).is_err());
        m.attachments.uncertainty = "widths".into();
        m.y = vec![1, 2];
        m.attachments.bands.insert(2, [Some(4), None]);
        assert!(import::prepare(&table, &m, || false).is_ok()); // deliberately shared supplied widths
    }
}
