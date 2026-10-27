use std::{
    path::Path,
    process::{Command, Output},
};
fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_1dPetePlot"))
        .args(args)
        .env_remove("DISPLAY")
        .env_remove("WAYLAND_DISPLAY")
        .output()
        .unwrap()
}
fn write(dir: &Path, name: &str, data: &str) -> String {
    let p = dir.join(name);
    std::fs::write(&p, data).unwrap();
    p.to_str().unwrap().into()
}
fn success(out: Output) {
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}
#[test]
fn real_inputs_export_all_four_modes_without_display() {
    let d = tempfile::tempdir().unwrap();
    let dir = d.path();
    let data = write(
        dir,
        "signals.1D",
        "# two series\n0.1 0.2\n0.4 0.3\n0.6 0.7\n0.3 0.4\n",
    );
    let x = write(dir, "x.1D", "0.2\n0.3\n0.5\n0.4\n");
    let matrix = write(dir, "rdm.1D", "0 0.2 0.8\n0.2 0 0.6\n0.8 0.6 0\n");
    let order = write(dir, "order.1D", "2\n0\n1\n");
    let output = dir.join("result.png");
    let output = output.to_str().unwrap();
    success(run(&["-infiles", &data, "-prefix", output]));
    success(run(&[
        "-histogram",
        "-histbin",
        "8",
        "-density",
        "-vline",
        "-0.1",
        "0",
        "-infiles",
        &format!("{data}[0]"),
        "-prefix",
        output,
        "-overwrite",
    ]));
    success(run(&[
        "-scatter",
        "-diagonal",
        "-xfile",
        &x,
        "-infiles",
        &format!("{data}[1]"),
        "-prefix",
        output,
        "-overwrite",
    ]));
    success(run(&[
        "-heat",
        "-zerocenter",
        "-cbar",
        "Reds_and_Blues_Inv",
        "-reorder",
        &order,
        "-infiles",
        &matrix,
        "-prefix",
        output,
        "-overwrite",
    ]));
    let png = std::fs::read(output).unwrap();
    assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    let svg = dir.join("result.svg");
    success(run(&[
        "-heat",
        "-infiles",
        &matrix,
        "-prefix",
        svg.to_str().unwrap(),
    ]));
    assert!(std::fs::read_to_string(svg).unwrap().contains("<svg"));
}
#[test]
fn validation_is_actionable_and_does_not_create_output() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("bad.png");
    let out = p.to_str().unwrap();
    for (args, message) in [
        (
            vec!["-demo", "-width", "900", "-dpi", "200"],
            "cannot be combined",
        ),
        (vec!["-demo", "-histbin", "0"], "histbin"),
        (vec!["-demo", "-vline_hl", "0.2"], "requires an argument"),
        (vec!["-demo", "-mystery"], "unsupported option"),
        (vec!["-demo", "-height", "1"], "figure size"),
        (vec!["-diagonal", "anything.1D"], "requires -scatter"),
    ] {
        let mut a = args;
        a.extend(["-prefix", out]);
        let result = run(&a);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains(message));
        assert!(!p.exists());
    }
}
#[test]
fn inches_dpi_and_overwrite_contract() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("figure.png");
    let path = p.to_str().unwrap();
    success(run(&[
        "-demo", "heat", "-figsize", "6", "5", "-dpi", "100", "-prefix", path,
    ]));
    let before = std::fs::read(&p).unwrap();
    assert_eq!(u32::from_be_bytes(before[16..20].try_into().unwrap()), 600);
    assert_eq!(u32::from_be_bytes(before[20..24].try_into().unwrap()), 500);
    let out = run(&["-demo", "-prefix", path]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("-overwrite"));
    assert_eq!(std::fs::read(&p).unwrap(), before);
}
#[test]
fn invalid_reorder_and_ragged_data_are_rejected() {
    let d = tempfile::tempdir().unwrap();
    let m = write(d.path(), "m.1D", "0 1\n1 0\n");
    let r = write(d.path(), "order.1D", "0\n0\n");
    let p = d.path().join("bad.svg");
    let out = run(&[
        "-heat",
        "-reorder",
        &r,
        "-infiles",
        &m,
        "-prefix",
        p.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("permutation"));
    assert!(!p.exists());
    write(d.path(), "m.1D", "1 2\n3\n");
    let out = run(&[&m, "-prefix", p.to_str().unwrap()]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("ragged"));
}
#[cfg(not(feature = "gui"))]
#[test]
fn headless_build_reports_how_to_export() {
    let out = run(&["-demo"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("no native-window support"));
}

#[test]
fn panel_title_and_arrangement_options_are_validated() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("arranged.svg");
    success(run(&[
        "-demo",
        "-title",
        "My RSA figure",
        "-panel-title",
        "1",
        "Evidence",
        "-panel-title",
        "4",
        "RDM",
        "-panel-order",
        "4,1,2,3",
        "-columns",
        "2",
        "-prefix",
        p.to_str().unwrap(),
    ]));
    for (args, message) in [
        (
            vec!["-demo", "-panel-order", "1,1,2,3"],
            "every panel exactly once",
        ),
        (
            vec!["-demo", "-panel-title", "5", "Missing"],
            "does not exist",
        ),
        (vec!["-demo", "-panel-title", "0", "Missing"], "start at 1"),
        (
            vec!["-demo", "-panel-title", "1", "A", "-panel-title", "1", "B"],
            "more than one title",
        ),
        (vec!["-demo", "-columns", "0"], "columns must be"),
        (vec!["-demo", "-columns", "5"], "columns must be"),
        (
            vec!["-demo", "-panel-order", "4,3,2,1", "-reverse_order"],
            "choose -panel-order",
        ),
    ] {
        let result = run(&args);
        assert!(!result.status.success());
        assert!(
            String::from_utf8_lossy(&result.stderr).contains(message),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
}

#[test]
fn apqc_patterns_export_bands_censor_boxplots_and_thin_strip() {
    let d = tempfile::tempdir().unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    let data = format!("{}[enorm,outlier]", source.join("qc.tsv").display());
    let pm = format!("{}[se_enorm,se_outlier]", source.join("qc.tsv").display());
    let censor = source.join("qc_censor.1D");
    let censor = censor.to_str().unwrap();
    for ext in ["png", "svg"] {
        let out = d.path().join(format!("qc.{ext}"));
        success(run(&[
            "-boxplot_on",
            "-patches",
            "60",
            "60",
            "-colors",
            "black",
            "-infiles",
            &data,
            "-yfiles_pm",
            &pm,
            "-ylim_use_pm",
            "-censor_files",
            censor,
            "-censor_trs",
            "0..2",
            "-censor_hline",
            "0.2",
            "0.05",
            "-yaxis",
            "0:0.6",
            "0:",
            "-ylabels",
            "enorm (~mm)",
            "outlier frac",
            "-title",
            "Synthetic QC",
            "-prefix",
            out.to_str().unwrap(),
        ]));
    }
    let out = d.path().join("strip.png");
    success(run(&[
        "-margin_off",
        "-one_graph",
        "-figsize",
        "12",
        "0.5",
        "-dpi",
        "150",
        "-patches",
        "60",
        "60",
        "-infiles",
        &data,
        "-scale",
        "SCALE_TO_HLINE",
        "-yaxis",
        "0:3",
        "-censor_files",
        censor,
        "-censor_hline",
        "0.2",
        "0.05",
        "-colors",
        "blue",
        "green",
        "-prefix",
        out.to_str().unwrap(),
    ]));
    let bytes = std::fs::read(out).unwrap();
    assert_eq!(u32::from_be_bytes(bytes[16..20].try_into().unwrap()), 1800);
    assert_eq!(u32::from_be_bytes(bytes[20..24].try_into().unwrap()), 75);
    let motion = source.join("qc_motion.1D");
    success(run(&[
        "-sepscl",
        "-boxplot_on",
        "-patches",
        "60",
        "60",
        "-reverse_order",
        "-infiles",
        motion.to_str().unwrap(),
        "-ylabels",
        "VOLREG",
        "-censor_files",
        censor,
        "-prefix",
        d.path().join("motion.png").to_str().unwrap(),
    ]));
}

#[test]
fn missing_data_exports_and_invalid_qc_leave_output_untouched() {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("gap.svg");
    let path = p.to_str().unwrap();
    for mode in ["gap", "mean"] {
        success(run(&[
            "-infiles",
            "1D: 1 2 NaN 4 3",
            "-missing",
            mode,
            "-prefix",
            path,
            "-overwrite",
        ]));
    }
    let before = std::fs::read(&p).unwrap();
    for args in [
        vec!["-infiles", "1D: 1 2 NaN 4 3"],
        vec!["-infiles", "1D: 1 2 3", "-censor_files", "1D: 1 0"],
        vec!["-infiles", "1D: 1 2 3", "-patches", "2", "2"],
        vec![
            "-infiles",
            "1D: 1 2 3",
            "-scale",
            "SCALE_TO_HLINE",
            "-censor_hline",
            "0",
        ],
        vec!["-infiles", "1D: 1 2 3", "-boxplot_on", "-one_graph"],
        vec!["-histogram", "-infiles", "1D: 1 2 3", "-censor_trs", "1"],
        vec!["-infiles", "1D: 1 2 3", "-yaxis", "5:0"],
    ] {
        let mut a = args;
        a.extend(["-prefix", path, "-overwrite"]);
        assert!(!run(&a).status.success());
        assert_eq!(std::fs::read(&p).unwrap(), before);
    }
}

#[test]
fn independent_and_matrix_shading_export_and_validate() {
    let d = tempfile::tempdir().unwrap();
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let data = root.join("examples/qc_motion.1D");
    let mat = root.join("examples/X.demo.mat.1D");
    for ext in ["png", "svg"] {
        success(run(&[
            "-infiles",
            &format!("{}[0]", data.display()),
            "-shade",
            "5:12",
            "-shade-color",
            "purple",
            "-shade-label",
            "Practice",
            "-shade-alpha",
            "0.3",
            "-shade-xmat",
            mat.to_str().unwrap(),
            "-shade-xmat-colors",
            "blue",
            "orange",
            "-boxplot_on",
            "-censor_trs",
            "0..2",
            "-prefix",
            d.path().join(format!("conditions.{ext}")).to_str().unwrap(),
        ]));
    }
    for args in [
        vec!["-shade-color", "blue"],
        vec!["-shade-label", "Oops"],
        vec!["-shade", "2:1"],
        vec!["-shade", "0:2", "-shade-alpha", "1.1"],
        vec!["-shade", "0:2", "-shade-color", "unknown"],
        vec!["-shade-xmat-threshold", "-0.1"],
        vec!["-shade-xmat-colors", "red"],
    ] {
        let p = d.path().join("bad.png");
        let mut a = args;
        a.extend(["-infiles", "1D: 1 2 3", "-prefix", p.to_str().unwrap()]);
        assert!(!run(&a).status.success());
        assert!(!d.path().join("bad.png").exists());
    }
}

#[test]
fn saved_project_replays_headlessly_with_caption_views_sizes_and_overwrite_policy() {
    let temp = tempfile::tempdir().unwrap();
    let recipe = temp.path().join("edited.peteplot.json");
    let mut project = serde_json::json!({
        "format": "1dPetePlot", "version": 1,
        "origin": {"arguments": ["-demo", "line"], "sources": []},
        "state": {
            "presentation": {"title": "Edited from the native window", "caption": "Caption with α and r^{2}.\nA second line survives replay.", "columns": 1,
                "panels": [{"id":0,"title":"Edited graph","caption":"Panel caption","xlabel":"Time (s)","ylabel":"Association"}]},
            "views": [{"x":[10,70],"y":[-0.1,0.8]}], "width": 1200, "height": 900, "full_plot": false
        }
    });
    std::fs::write(&recipe, serde_json::to_vec_pretty(&project).unwrap()).unwrap();
    let png = temp.path().join("replay.png");
    let svg = temp.path().join("replay.svg");
    let recipe = recipe.to_str().unwrap();
    let png = png.to_str().unwrap();
    let svg = svg.to_str().unwrap();
    success(run(&["-recipe", recipe, "-prefix", png]));
    let original = std::fs::read(png).unwrap();
    assert_eq!(
        u32::from_be_bytes(original[16..20].try_into().unwrap()),
        1200
    );
    assert_eq!(
        u32::from_be_bytes(original[20..24].try_into().unwrap()),
        900
    );
    let denied = run(&["-recipe", recipe, "-prefix", png]);
    assert!(!denied.status.success());
    assert_eq!(std::fs::read(png).unwrap(), original);
    success(run(&["-recipe", recipe, "-prefix", png, "-overwrite"]));
    assert_eq!(std::fs::read(png).unwrap(), original);
    success(run(&["-recipe", recipe, "-prefix", svg, "-width", "1600"]));
    assert!(
        std::fs::read_to_string(svg)
            .unwrap()
            .contains("width=\"1600\" height=\"900\"")
    );
    assert!(
        !run(&[
            "-recipe",
            recipe,
            "-title",
            "Unrecorded override",
            "-prefix",
            png,
            "-overwrite"
        ])
        .status
        .success()
    );
    project["version"] = serde_json::json!(99);
    std::fs::write(recipe, serde_json::to_vec(&project).unwrap()).unwrap();
    assert!(
        !run(&["-recipe", recipe, "-prefix", png, "-overwrite"])
            .status
            .success()
    );
    assert_eq!(std::fs::read(png).unwrap(), original);
}

#[test]
fn mixed_version_two_recipe_replays_imports_and_rejects_invalid_identities() {
    let temp = tempfile::tempdir().unwrap();
    let recipe = temp.path().join("mixed.peteplot.json");
    let mut project = serde_json::json!({
        "format": "1dPetePlot", "version": 2,
        "origin": {"arguments": [], "sources": [], "imports": [
            {"id": 4, "origin": {"arguments": ["-infiles", "1D: 1 3 2 4", "-one"], "sources": []},
             "roles": {"columns": ["effect"], "x": null, "y": [0], "kind": "Lines", "missing": "error"}},
            {"id": 9, "origin": {"arguments": ["-infiles", "1D: 1 1 2 3", "-histogram", "-one"], "sources": []},
             "roles": {"columns": ["effect"], "x": null, "y": [0], "kind": "Distribution", "missing": "error"}}
        ]},
        "state": {
            "presentation": {"title": "Mixed imports", "subtitle": "Independent sources", "caption": "Saved presentation", "columns": 2,
                "panels": [
                    {"id":9,"title":"Distribution","caption":"Full counts","xlabel":"Effect","ylabel":"Count"},
                    {"id":4,"title":"Curve","caption":"Original observations","xlabel":"Row","ylabel":"Effect"}]},
            "views": [{"x":[0,4],"y":[0,3]}, {"x":[0,3],"y":[0,5]}], "width": 1200, "height": 700, "full_plot": true
        }
    });
    std::fs::write(&recipe, serde_json::to_vec(&project).unwrap()).unwrap();
    for ext in ["png", "svg"] {
        success(run(&[
            "-recipe",
            recipe.to_str().unwrap(),
            "-prefix",
            temp.path().join(format!("mixed.{ext}")).to_str().unwrap(),
        ]));
    }
    let output = temp.path().join("bad.png");
    let valid = project.clone();
    for bad in ["identity", "roles", "nested"] {
        project = valid.clone();
        match bad {
            "identity" => project["origin"]["imports"][1]["id"] = serde_json::json!(4),
            "roles" => project["origin"]["imports"][1]["roles"]["y"] = serde_json::json!([2]),
            _ => {
                project["origin"]["imports"][0]["origin"]["imports"] =
                    valid["origin"]["imports"].clone()
            }
        }
        std::fs::write(&recipe, serde_json::to_vec(&project).unwrap()).unwrap();
        assert!(
            !run(&[
                "-recipe",
                recipe.to_str().unwrap(),
                "-prefix",
                output.to_str().unwrap()
            ])
            .status
            .success()
        );
        assert!(!output.exists());
    }
    project = valid;
    project["origin"]["imports"] = serde_json::json!([]);
    project["state"]["presentation"]["panels"] = serde_json::json!([]);
    project["state"]["presentation"]["columns"] = serde_json::json!(1);
    project["state"]["views"] = serde_json::json!([]);
    std::fs::write(&recipe, serde_json::to_vec(&project).unwrap()).unwrap();
    let out = run(&[
        "-recipe",
        recipe.to_str().unwrap(),
        "-prefix",
        output.to_str().unwrap(),
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("add a graph"));
    assert!(!output.exists());
}

#[test]
fn rsa_presets_export_without_display_and_invalid_schemas_preserve_outputs() {
    let d = tempfile::tempdir().unwrap();
    let source = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/rsa/ordinary.rsa.1D");
    let file = source.to_str().unwrap();
    let png = d.path().join("rsa.png");
    let svg = d.path().join("rsa.svg");
    success(run(&["-rsa", file, "-prefix", png.to_str().unwrap()]));
    success(run(&[
        "-rsa",
        file,
        "-rsa-rois",
        "12,3,27",
        "-rsa-models",
        "visual,semantic",
        "-rsa-plots",
        "scatter,heatmap",
        "-rsa-pair",
        "semantic,visual",
        "-prefix",
        svg.to_str().unwrap(),
    ]));
    assert!(std::fs::read_to_string(svg).unwrap().contains("<svg"));
    let before = std::fs::read(&png).unwrap();
    for args in [
        vec!["-rsa-effect", "pr"],
        vec!["-rsa-rois", "9999"],
        vec!["-rsa-models", "fit"],
        vec!["-rsa-version", "2"],
        vec!["-missing", "mean"],
        vec!["-xfile", file],
        vec!["-rsa-plots", "scatter", "-rsa-models", "semantic"],
        vec!["-rsa-plots", "heatmap", "-histbin", "5"],
    ] {
        let mut cmd = vec!["-rsa", file];
        cmd.extend(args);
        cmd.extend(["-prefix", png.to_str().unwrap(), "-overwrite"]);
        assert!(!run(&cmd).status.success());
        assert_eq!(std::fs::read(&png).unwrap(), before);
    }
    let regression = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/rsa/regression.rsa.1D");
    success(run(&[
        "-rsa",
        regression.to_str().unwrap(),
        "-rsa-effect",
        "pr",
        "-prefix",
        png.to_str().unwrap(),
        "-overwrite",
    ]));
}
