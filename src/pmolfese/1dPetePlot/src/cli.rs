use anyhow::{Context, Result, bail, ensure};

pub const HELP: &str = r#"1dPetePlot — native & headless plotting prototype

Usage: 1dPetePlot [options] data.1D ...
       1dPetePlot -demo [line|histogram|scatter|heat|gallery] [options]
       1dPetePlot                    empty native workspace (GUI build)

Window:    Open Data (Cmd/Ctrl-I) or drop .1D/headered TSV files; preview, then Add.
           Choose numeric columns and missing-value policy before adding a graph.

Input:     -infiles/-yfiles FILE ...   (or positional files; '-' / -stdin for stdin)
           -xfile/-x FILE             explicit finite x values
           [cols]{rows} selectors; TSV [column_name]; trailing ' transposes
           '1D: 3@0,2@1 | 1 2 3 4 5' inline series (quote in your shell)
           -missing error|gap|mean    strict default; gaps or column-mean imputation
Plots:     -scatter [-diagonal]       points; identity line with equal scaling
           -histogram [-histbin N | -histwidth W] [-density]
           -heat [-reorder FILE] [-zerocenter] [-cbar NAME]
           -one_graph/-one | -sep     overlay or separate panels (default separate)
           -sepscl                    independent y ranges
RSA:       -rsa RESULTS.rsa.1D       ROI-wide, time-series or run results; preset v1
           -rsa-family effects|fit|weights|fit-contrast|commonality|partial-r2
           -rsa-measures COL,COL    explicit named quantities in a non-effects family
           -rsa-models A,B           select up to 8 models (default all)
           -rsa-effect auto|r|b|pr|effect   auto: r (wide), effect (long), b (regression)
           -rsa-rois 3,12,27         ROI IDs in chosen order (default all, ID order)
           -rsa-plots histogram,scatter,heatmap,intervals
           Long tables: curves/intervals (time), points (runs), heatmap (both)
           -rsa-times 0,2,5          Select time_index IDs; retain index order
           -rsa-summaries MEAN,A,B   Select run/contrast summaries in chosen order
           -rsa-ci auto|boot|cboot|dual   Supplied CI method for intervals (auto requires one).
                                        Up to 40 ROIs; regression CIs are beta, never partial r.
           -rsa-pair XMODEL,YMODEL   one comparison (default first two models)
           Wide default: distributions, one comparison, heatmap if <=40 ROIs.
           Long default: curves (time), categorical points (runs), one panel per ROI.
           Explicit -rsa-ci boot|cboot|dual also adds supplied bands to time curves.
           -missing gap omits missing effects/pairs; RSA imputation is rejected.
           For intervals, incomplete effect/CI triplets retain a row marked missing.
           CI levels come from source metadata; missing levels are never assumed.
           Recognized RSA files also work with Open Data / drop-to-preview.
RDM:       -rdm MATRIX.1D ...         1..8 saved square symmetric matrices
           -rdm-labels FILE          one unique condition name per line in source order;
                                     one shared file, or repeat once per matrix
           -rdm-order FILE           desired condition names, one per line (both axes)
           -rdm-kind unspecified|similarity|dissimilarity (default unspecified)
           -rdm-scale independent|shared (default independent; shared is opt-in)
           -clim MIN:MAX             fixed common color limits; -cbar/-zerocenter supported
           Labels must match exactly; negatives/diagonal preserved, no conversion.
           Open/drop a numeric .1D matrix, choose RDM / named matrix, attach labels.
Project:   -recipe FILE.peteplot.json [-prefix FILE.png|FILE.svg] [-show]
           Replay saved edits, panel order, views, and dimensions.
Drawing:   -caption TEXT -panel-caption N TEXT (repeatable; 1-based)
           -title/-plabel TEXT -xlabel TEXT -ylabels TEXT ...
           -panel-title N TEXT        title of original panel N (repeatable; 1-based)
           -panel-order 3,1,2         all panel numbers, in display order
           -columns N                 panels per row (1..4)
           -legend_labels TEXT ... -legend_on -reverse_order
           -xaxis MIN:MAX -yaxis MIN:MAX ... (also MIN:, :MAX, :)
           -vline X ... -vline_hl X LABEL (repeatable)
           -colors '#RRGGBB' ...      names, #RGB, grayscale 0..1
           -ylabels VOLREG|ALLINPAR6|ALLINPAR9|ALLINPAR12 ...
           -ylabels_maxlen N          wrap axis labels
           Labels: Unicode Greek; x_{i} and r^{2}; bare underscores stay literal
Shading:   -shade START:END           x-coordinate interval on all panels; repeatable
           -shade-color COLOR        sets the most recent range (default palette)
           -shade-label TEXT         labeled color key; does not censor samples
           -shade-alpha A            opacity 0..1 (default 0.18)
           -shade-xmat X.mat.1D       positive condition groups; full plotted timeline
           -shade-xmat-colors COLOR ...   condition palette override
           -shade-xmat-threshold F    fraction of condition peak |response| (default .05)
           Matrix shading shows modeled activity, NOT exact stimulus timing.
Grouping:  -group-by COLUMN -one      original column name/index; no aggregation
           groups retain original rows; at most 16 series; no mean imputation
QC/lines:  -yfiles_lo FILE ... -yfiles_hi FILE ...   supplied asymmetric endpoints
           -yfiles_pm FILE ... [-ylim_use_pm]  symmetric bands, absolute widths
           -censor_files FILE ...     0=censored, 1=retained; -censor FILE alias
           -censor_trs '0..3,8,$' ...  union of selected-row indices and masks
           -censor_RGB 'red 0.3'      color [opacity] or 'R G B [A]'
           -censor_hline VALUE|NONE ...  one threshold or one per series
           -patches N1 N2 ...         run lengths summing to selected row count
           -boxplot_on [-bplot_view BC_ONLY|AC_ONLY]
           -scale VALUE ...          or SCALE_TO_HLINE / SCALE_TO_MAX
           -margin_off               one-panel plot-only strip (no boxplots)
           -dx/-dt STEP -xzero VALUE
           -cbar viridis|Reds_and_Blues|Reds_and_Blues_Inv
           -clim MIN:MAX             explicit heatmap color limits
Output:    -prefix FILE.png|FILE.svg  save and exit; extension defaults to .png
           -png FILE / -svg FILE     output aliases
           -width W -height H        exact pixels (default 1500x1200, grows for more panels)
           -figsize W H -dpi D        inches and DPI (default 150); not with pixels
           -overwrite                permit replacing a CLI output
Window:    -show                     also show after saving
           drag to pan, Shift-drag a zoom box, double-click or R to reset
           right-click to copy one graph or all; Save opens a native file dialog
           Save exports the current view; full plot is selectable
Other:     -help/-h/--help -version

No output filename opens a native window. PNG/SVG export needs no window or GPU.
Native macOS / Linux Wayland only; no X11 or XQuartz. Build --no-default-features
for a CLI-only executable. This prototype does not implement every legacy flag;
unknown options fail explicitly. Missing values require an explicit policy.
Demo data are synthetic, not experimental results. PDF and full legacy parity
are later milestones. See README.md for supported syntax and limitations.
"#;

#[derive(Clone, Debug)]
pub struct Options {
    pub arguments: Vec<String>,
    pub source_args: Vec<usize>,
    pub rsa: crate::rsa::Request,
    pub rdm: crate::rdm::Request,
    pub caption: String,
    pub panel_captions: Vec<(usize, String)>,
    pub shades: Vec<crate::shading::Shade>,
    pub shade_design: crate::shading::Design,
    pub qc: crate::qc::Options,
    pub missing: crate::input::Missing,
    pub auto_size: bool,
    pub files: Vec<String>,
    pub xfile: Option<String>,
    pub group_by: Option<String>,
    pub demo: Option<String>,
    pub kind: String,
    pub one: bool,
    pub sepscl: bool,
    pub reverse: bool,
    pub output: Option<String>,
    pub show: bool,
    pub overwrite: bool,
    pub width: u32,
    pub height: u32,
    pub title: Option<String>,
    pub panel_titles: Vec<(usize, String)>,
    pub panel_order: Option<Vec<usize>>,
    pub columns: Option<usize>,
    pub xlabel: String,
    pub ylabels: Vec<String>,
    pub legends: Vec<String>,
    pub legend_on: bool,
    pub vlines: Vec<(f64, Option<String>)>,
    pub bins: Option<usize>,
    pub binwidth: Option<f64>,
    pub density: bool,
    pub diagonal: bool,
    pub zerocenter: bool,
    pub cmap: String,
    pub reorder: Option<String>,
    pub clim: Option<[f64; 2]>,
    pub xlim: Option<[f64; 2]>,
    pub ylim: Option<[f64; 2]>,
    pub colors: Vec<String>,
    pub dx: f64,
    pub xzero: f64,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            arguments: vec![],
            source_args: vec![],
            rsa: Default::default(),
            rdm: Default::default(),
            caption: String::new(),
            panel_captions: vec![],
            shades: vec![],
            shade_design: Default::default(),
            qc: Default::default(),
            missing: Default::default(),
            auto_size: true,
            files: vec![],
            xfile: None,
            group_by: None,
            demo: None,
            kind: "line".into(),
            one: false,
            sepscl: false,
            reverse: false,
            output: None,
            show: false,
            overwrite: false,
            width: 1500,
            height: 1200,
            title: None,
            panel_titles: vec![],
            panel_order: None,
            columns: None,
            xlabel: String::new(),
            ylabels: vec![],
            legends: vec![],
            legend_on: false,
            vlines: vec![],
            bins: None,
            binwidth: None,
            density: false,
            diagonal: false,
            zerocenter: false,
            cmap: "viridis".into(),
            reorder: None,
            clim: None,
            xlim: None,
            ylim: None,
            colors: vec![],
            dx: 1.,
            xzero: 0.,
        }
    }
}
fn number(s: &str) -> Result<f64> {
    let n: f64 = s
        .parse()
        .with_context(|| format!("expected a number, got '{s}'"))?;
    ensure!(n.is_finite(), "expected a finite number, got '{s}'");
    Ok(n)
}
fn limits(s: &str) -> Result<[f64; 2]> {
    let (a, b) = s.split_once(':').context("range must be MIN:MAX")?;
    let v = [number(a)?, number(b)?];
    ensure!(v[0] < v[1], "range minimum must be less than maximum");
    Ok(v)
}
fn is_numeric_arg(s: &str) -> bool {
    s.parse::<f64>().is_ok()
        || s.split_once(':').is_some_and(|(a, b)| {
            a.parse::<f64>().is_ok() && (b.is_empty() || b.parse::<f64>().is_ok())
        })
}
fn argument<'a>(args: &'a [String], i: &mut usize) -> Result<&'a str> {
    let flag = &args[*i];
    *i += 1;
    let value = args
        .get(*i)
        .map(String::as_str)
        .with_context(|| format!("{flag} requires an argument"))?;
    ensure!(
        !value.starts_with('-') || value == "-" || is_numeric_arg(value),
        "{flag} requires an argument before option '{value}'"
    );
    Ok(value)
}
fn list(args: &[String], i: &mut usize) -> Result<Vec<String>> {
    let flag = args[*i].clone();
    let mut out = vec![];
    while let Some(v) = args.get(*i + 1) {
        if v.starts_with('-') && v != "-" && !is_numeric_arg(v) {
            break;
        }
        *i += 1;
        out.push(v.clone());
    }
    ensure!(!out.is_empty(), "{flag} requires at least one value");
    Ok(out)
}
pub fn parse(args: &[String]) -> Result<Options> {
    let mut o = Options {
        arguments: args.to_vec(),
        ..Default::default()
    };
    let mut i = 0;
    let (mut pixel, mut inches, mut dpi_used) = (false, None, false);
    let mut dpi = 150.;
    let mut kind_set = false;
    let mut used_flags = vec![];
    while i < args.len() {
        let start = i;
        let flag = args[i].as_str();
        used_flags.push(flag);
        match flag {
            "-rdm" => o.rdm.files.extend(list(args, &mut i)?),
            "-rdm-labels" => o.rdm.labels.push(argument(args, &mut i)?.into()),
            "-rdm-order" => o.rdm.order = Some(argument(args, &mut i)?.into()),
            "-rdm-scale" => o.rdm.scale = argument(args, &mut i)?.into(),
            "-rdm-kind" => o.rdm.sense = argument(args, &mut i)?.into(),
            "-rdm-version" => o.rdm.version = argument(args, &mut i)?.parse()?,
            "-rsa" => {
                ensure!(o.rsa.file.is_none(), "specify one RSA table per preset");
                o.rsa.file = Some(argument(args, &mut i)?.into());
            }
            "-rsa-family" => o.rsa.family = argument(args, &mut i)?.into(),
            "-rsa-measures" => o.rsa.measures = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-models" => o.rsa.models = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-rois" => {
                o.rsa.rois = Some(
                    crate::rsa::names(argument(args, &mut i)?)?
                        .iter()
                        .map(|v| v.parse())
                        .collect::<std::result::Result<_, _>>()?,
                )
            }
            "-rsa-plots" => o.rsa.plots = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-effect" => o.rsa.effect = argument(args, &mut i)?.into(),
            "-rsa-times" => o.rsa.times = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-summaries" => o.rsa.summaries = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-ci" => o.rsa.ci = argument(args, &mut i)?.into(),
            "-rsa-pair" => o.rsa.pair = Some(crate::rsa::names(argument(args, &mut i)?)?),
            "-rsa-version" => o.rsa.version = argument(args, &mut i)?.parse()?,
            "-caption" => o.caption = argument(args, &mut i)?.into(),
            "-panel-caption" => {
                let n = argument(args, &mut i)?.parse::<usize>()?;
                o.panel_captions.push((n, argument(args, &mut i)?.into()));
            }
            "-infiles" | "-yfiles" => o.files.extend(list(args, &mut i)?),
            "-stdin" => o.files.push("-".into()),
            "-xfile" | "-x" => o.xfile = Some(argument(args, &mut i)?.into()),
            "-demo" => {
                o.demo = Some("gallery".into());
                if let Some(v) = args.get(i + 1)
                    && !v.starts_with('-')
                {
                    i += 1;
                    o.demo = Some(v.clone());
                }
            }
            "-scatter" | "-heat" | "-histogram" => {
                ensure!(!kind_set, "choose only one plot type");
                kind_set = true;
                o.kind = flag[1..].into();
            }
            "-prefix" | "-png" | "-svg" => {
                ensure!(o.output.is_none(), "specify only one output filename");
                let mut name = argument(args, &mut i)?.to_owned();
                if std::path::Path::new(&name).extension().is_none() {
                    name.push_str(if flag == "-svg" { ".svg" } else { ".png" });
                }
                if flag == "-png" || flag == "-svg" {
                    ensure!(
                        name.to_lowercase().ends_with(&format!(".{}", &flag[1..])),
                        "output extension conflicts with {flag}"
                    );
                }
                o.output = Some(name);
            }
            "-show" => o.show = true,
            "-overwrite" => o.overwrite = true,
            "-one" | "-one_graph" => o.one = true,
            "-sep" => o.one = false,
            "-sepscl" => o.sepscl = true,
            "-reverse_order" => o.reverse = true,
            "-width" | "-height" => {
                pixel = true;
                let v: u32 = argument(args, &mut i)?
                    .parse()
                    .context("size must be a positive integer")?;
                if flag == "-width" {
                    o.width = v
                } else {
                    o.height = v
                }
            }
            "-figsize" => {
                let w = number(argument(args, &mut i)?)?;
                let h = number(argument(args, &mut i)?)?;
                ensure!(w > 0. && h > 0., "figure dimensions must be positive");
                inches = Some((w, h));
            }
            "-dpi" => {
                dpi_used = true;
                dpi = number(argument(args, &mut i)?)?;
                ensure!(dpi > 0., "DPI must be positive");
            }
            "-title" | "-plabel" => o.title = Some(argument(args, &mut i)?.into()),
            "-panel-title" => {
                let n: usize = argument(args, &mut i)?
                    .parse()
                    .context("panel number must be a positive integer")?;
                ensure!(n > 0, "panel numbers start at 1");
                let title = argument(args, &mut i)?.into();
                ensure!(
                    !o.panel_titles.iter().any(|(i, _)| *i == n),
                    "panel {n} has more than one title"
                );
                o.panel_titles.push((n, title));
            }
            "-panel-order" => {
                ensure!(o.panel_order.is_none(), "specify -panel-order only once");
                o.panel_order = Some(
                    argument(args, &mut i)?
                        .split(',')
                        .map(|v| v.trim().parse::<usize>())
                        .collect::<std::result::Result<Vec<_>, _>>()
                        .context("panel order must be comma-separated numbers, e.g. 3,1,2")?,
                );
            }
            "-columns" => {
                let n: usize = argument(args, &mut i)?
                    .parse()
                    .context("columns must be 1..4")?;
                ensure!((1..=4).contains(&n), "columns must be 1..4");
                o.columns = Some(n);
            }
            "-xlabel" => o.xlabel = argument(args, &mut i)?.into(),
            "-ylabel" | "-ylabels" => {
                o.ylabels = list(args, &mut i)?
                    .into_iter()
                    .flat_map(|v| {
                        if v == "VOLREG" {
                            vec![
                                "roll (deg)",
                                "pitch (deg)",
                                "yaw (deg)",
                                "dS (mm)",
                                "dL (mm)",
                                "dP (mm)",
                            ]
                            .into_iter()
                            .map(String::from)
                            .collect()
                        } else if ["ALLINPAR6", "ALLINPAR9", "ALLINPAR12"].contains(&v.as_str()) {
                            let n = v[8..].parse::<usize>().unwrap();
                            [
                                "x-shift (mm)",
                                "y-shift (mm)",
                                "z-shift (mm)",
                                "z-angle (deg)",
                                "x-angle (deg)",
                                "y-angle (deg)",
                                "x-scale",
                                "y-scale",
                                "z-scale",
                                "y/x-shear",
                                "z/x-shear",
                                "z/y-shear",
                            ][..n]
                                .iter()
                                .map(|s| s.to_string())
                                .collect()
                        } else {
                            vec![v]
                        }
                    })
                    .collect();
            }
            "-legend_labels" => {
                o.legends = list(args, &mut i)?;
                o.legend_on = true;
            }
            "-legend_on" => o.legend_on = true,
            "-vline" => {
                for v in list(args, &mut i)? {
                    o.vlines.push((number(&v)?, None));
                }
            }
            "-vline_hl" => {
                let x = number(argument(args, &mut i)?)?;
                let label = argument(args, &mut i)?.into();
                o.vlines.push((x, Some(label)));
            }
            "-histbin" => {
                let n = argument(args, &mut i)?.parse()?;
                ensure!((1..=10000).contains(&n), "histbin must be 1..10000");
                o.bins = Some(n);
            }
            "-histwidth" => {
                let w = number(argument(args, &mut i)?)?;
                ensure!(w > 0., "histwidth must be positive");
                o.binwidth = Some(w);
            }
            "-density" => o.density = true,
            "-diagonal" => o.diagonal = true,
            "-zerocenter" => o.zerocenter = true,
            "-cbar" => o.cmap = argument(args, &mut i)?.into(),
            "-reorder" => o.reorder = Some(argument(args, &mut i)?.into()),
            "-volreg" => {
                o.ylabels = [
                    "roll (deg)",
                    "pitch (deg)",
                    "yaw (deg)",
                    "dS (mm)",
                    "dL (mm)",
                    "dP (mm)",
                ]
                .map(String::from)
                .to_vec()
            }
            "-xaxis" => o.xlim = Some(limits(argument(args, &mut i)?)?),
            "-yaxis" => {
                o.qc.ylimits = list(args, &mut i)?
                    .iter()
                    .map(|s| {
                        let (a, b) = s
                            .split_once(':')
                            .context("yaxis range must be MIN:MAX, MIN:, :MAX, or :")?;
                        Ok([
                            if a.is_empty() { None } else { Some(number(a)?) },
                            if b.is_empty() { None } else { Some(number(b)?) },
                        ])
                    })
                    .collect::<Result<_>>()?;
                o.ylim = if o.qc.ylimits.len() == 1 {
                    match o.qc.ylimits[0] {
                        [Some(a), Some(b)] => Some([a, b]),
                        _ => None,
                    }
                } else {
                    None
                };
            }
            "-shade" => {
                let range = limits(argument(args, &mut i)?)?;
                o.shades
                    .push(crate::shading::Shade::new(range, o.shades.len()));
            }
            "-shade-color" => {
                let color = crate::model::color(argument(args, &mut i)?)?;
                o.shades
                    .last_mut()
                    .context("-shade-color must follow a -shade range")?
                    .color = color;
            }
            "-shade-label" => {
                let label = argument(args, &mut i)?.to_string();
                ensure!(!label.trim().is_empty(), "shade label must not be empty");
                o.shades
                    .last_mut()
                    .context("-shade-label must follow a -shade range")?
                    .label = Some(label);
            }
            "-shade-alpha" => {
                let a = number(argument(args, &mut i)?)?;
                ensure!((0.0..=1.0).contains(&a), "shade alpha must be 0..1");
                o.shades
                    .last_mut()
                    .context("-shade-alpha must follow a -shade range")?
                    .alpha = a;
            }
            "-shade-xmat" => {
                ensure!(
                    o.shade_design.file.is_none(),
                    "specify shade-xmat only once"
                );
                o.shade_design.file = Some(argument(args, &mut i)?.into());
            }
            "-shade-xmat-colors" => o.shade_design.colors = list(args, &mut i)?,
            "-shade-xmat-threshold" => {
                let t = number(argument(args, &mut i)?)?;
                ensure!(
                    (0.0..1.0).contains(&t),
                    "shade-xmat-threshold must be in [0,1)"
                );
                o.shade_design.threshold = t;
            }
            "-group-by" => o.group_by = Some(argument(args, &mut i)?.into()),
            "-yfiles_lo" => o.qc.lower.extend(list(args, &mut i)?),
            "-yfiles_hi" => o.qc.upper.extend(list(args, &mut i)?),
            "-yfiles_pm" => o.qc.pm.extend(list(args, &mut i)?),
            "-ylim_use_pm" => o.qc.use_pm = true,
            "-censor_files" => o.qc.censor_files.extend(list(args, &mut i)?),
            "-censor" => o.qc.censor_files.push(argument(args, &mut i)?.into()),
            "-censor_trs" => o.qc.censor_trs.extend(list(args, &mut i)?),
            "-censor_RGB" => o.qc.censor_color = list(args, &mut i)?,
            "-censor_hline" => {
                o.qc.hlines = list(args, &mut i)?
                    .iter()
                    .map(|s| {
                        if s == "NONE" {
                            Ok(None)
                        } else {
                            number(s).map(Some)
                        }
                    })
                    .collect::<Result<_>>()?
            }
            "-patches" => {
                o.qc.patches = list(args, &mut i)?
                    .iter()
                    .map(|s| s.parse())
                    .collect::<std::result::Result<_, _>>()?
            }
            "-boxplot_on" => o.qc.boxplot = true,
            "-bplot_view" => {
                let v = argument(args, &mut i)?;
                ensure!(
                    ["BC_ONLY", "AC_ONLY"].contains(&v),
                    "bplot_view must be BC_ONLY or AC_ONLY"
                );
                o.qc.box_view = v.into();
            }
            "-scale" => o.qc.scales = list(args, &mut i)?,
            "-margin_off" => o.qc.margin_off = true,
            "-ylabels_maxlen" => {
                let n = argument(args, &mut i)?.parse::<usize>()?;
                ensure!(n > 0, "ylabels_maxlen must be positive");
                o.qc.label_wrap = Some(n);
            }
            "-missing" => {
                o.missing = match argument(args, &mut i)? {
                    "error" => crate::input::Missing::Error,
                    "gap" => crate::input::Missing::Gap,
                    "mean" => crate::input::Missing::Mean,
                    _ => bail!("missing policy must be error, gap, or mean"),
                }
            }
            "-clim" => o.clim = Some(limits(argument(args, &mut i)?)?),
            "-colors" => o.colors = list(args, &mut i)?,
            "-dx" | "-dt" => o.dx = number(argument(args, &mut i)?)?,
            "-xzero" => o.xzero = number(argument(args, &mut i)?)?,
            "--" => {
                o.source_args.extend(i + 1..args.len());
                o.files.extend_from_slice(&args[i + 1..]);
                break;
            }
            "-" => o.files.push("-".into()),
            _ if flag.starts_with('-') => {
                bail!("unsupported option '{flag}'; see -help for prototype coverage")
            }
            _ => o.files.push(flag.into()),
        }
        match flag {
            "-infiles" | "-yfiles" | "-yfiles_pm" | "-yfiles_lo" | "-yfiles_hi"
            | "-censor_files" | "-xfile" | "-x" | "-reorder" | "-censor" | "-shade-xmat"
            | "-rsa" | "-rdm" | "-rdm-labels" | "-rdm-order" => o.source_args.extend(start + 1..=i),
            "-stdin" | "-" => o.source_args.push(start),
            _ if !flag.starts_with('-') => o.source_args.push(start),
            _ => {}
        }
        i += 1;
    }
    if !o.rdm.files.is_empty() {
        ensure!(
            o.files.is_empty(),
            "use -rdm FILE ... for all matrix inputs"
        );
        for flag in &used_flags {
            ensure!(
                [
                    "-rdm",
                    "-rdm-labels",
                    "-rdm-order",
                    "-rdm-scale",
                    "-rdm-kind",
                    "-rdm-version",
                    "-prefix",
                    "-png",
                    "-svg",
                    "-show",
                    "-overwrite",
                    "-width",
                    "-height",
                    "-figsize",
                    "-dpi",
                    "-title",
                    "-plabel",
                    "-caption",
                    "-panel-title",
                    "-panel-caption",
                    "-panel-order",
                    "-columns",
                    "-reverse_order",
                    "-missing",
                    "-cbar",
                    "-clim",
                    "-zerocenter"
                ]
                .contains(flag),
                "{flag} is not supported with -rdm; condition order uses -rdm-order"
            );
        }
        ensure!(
            o.missing != crate::input::Missing::Mean,
            "RDMs do not impute values"
        );
    } else {
        ensure!(
            !used_flags.iter().any(|f| f.starts_with("-rdm-")),
            "RDM options require -rdm FILE ..."
        );
    }
    if o.rsa.file.is_some() {
        for flag in &used_flags {
            ensure!(
                [
                    "-rsa",
                    "-rsa-models",
                    "-rsa-family",
                    "-rsa-measures",
                    "-rsa-rois",
                    "-rsa-plots",
                    "-rsa-ci",
                    "-rsa-times",
                    "-rsa-summaries",
                    "-rsa-effect",
                    "-rsa-pair",
                    "-rsa-version",
                    "-prefix",
                    "-png",
                    "-svg",
                    "-show",
                    "-overwrite",
                    "-width",
                    "-height",
                    "-figsize",
                    "-dpi",
                    "-title",
                    "-plabel",
                    "-caption",
                    "-panel-title",
                    "-panel-caption",
                    "-panel-order",
                    "-columns",
                    "-reverse_order",
                    "-histbin",
                    "-histwidth",
                    "-density",
                    "-missing",
                    "-colors",
                ]
                .contains(flag),
                "{flag} is not supported with -rsa; use -rsa-plots and -rsa-effect for scientific roles"
            );
        }
        ensure!(
            o.missing != crate::input::Missing::Mean,
            "RSA presets do not impute effects; use -missing gap to omit missing observations explicitly"
        );
        ensure!(
            o.rsa.version == 1,
            "unsupported RSA preset version (supported: 1)"
        );
    } else {
        ensure!(
            !used_flags.iter().any(|f| f.starts_with("-rsa-")),
            "RSA selection options require -rsa FILE"
        );
    }
    o.auto_size = !pixel && inches.is_none() && !dpi_used;
    ensure!(
        !(pixel && (inches.is_some() || dpi_used)),
        "pixel sizes cannot be combined with -figsize or -dpi"
    );
    ensure!(
        !(o.reverse && o.panel_order.is_some()),
        "choose -panel-order or -reverse_order"
    );
    if inches.is_some() || dpi_used {
        let (w, h) = inches.unwrap_or((10., 8.));
        o.width = (w * dpi).round() as u32;
        o.height = (h * dpi).round() as u32;
    }
    ensure!(
        (400..=8192).contains(&o.width)
            && (if o.qc.margin_off { 32 } else { 300 }..=8192).contains(&o.height),
        "figure size must be 400..8192 wide and 300..8192 high"
    );
    ensure!(
        o.width as u64 * o.height as u64 <= 32_000_000,
        "figure exceeds 32 million pixels"
    );
    ensure!(
        !(o.bins.is_some() && o.binwidth.is_some()),
        "choose -histbin or -histwidth"
    );
    ensure!(
        !o.diagonal || o.kind == "scatter",
        "-diagonal requires -scatter"
    );
    ensure!(
        o.rsa.file.is_some()
            || o.kind == "histogram"
            || (o.bins.is_none() && o.binwidth.is_none() && !o.density),
        "histogram options require -histogram"
    );
    ensure!(
        !o.rdm.files.is_empty()
            || o.kind == "heat"
            || (o.reorder.is_none() && !o.zerocenter && o.clim.is_none() && o.cmap == "viridis"),
        "heatmap options require -heat"
    );
    ensure!(
        o.kind != "heat" || (o.xfile.is_none() && !o.one && o.vlines.is_empty()),
        "heatmaps cannot use xfile, overlay, or vertical reference lines"
    );
    ensure!(
        o.kind != "histogram" || o.xfile.is_none(),
        "histogram bins supply x values; do not use -xfile"
    );
    ensure!(
        !(o.clim.is_some() && o.zerocenter),
        "choose explicit -clim or -zerocenter"
    );
    ensure!(
        o.demo.is_none() || (o.files.is_empty() && o.xfile.is_none()),
        "demo cannot be combined with input files"
    );
    ensure!(
        !o.rdm.files.is_empty() || o.rsa.file.is_some() || o.demo.is_some() || !o.files.is_empty(),
        "supply a .1D file or -demo (see -help)"
    );
    ensure!(
        o.group_by.is_none()
            || (o.demo.is_none()
                && o.files.len() == 1
                && o.one
                && o.kind != "heat"
                && o.rsa.file.is_none()
                && o.rdm.files.is_empty()),
        "group-by requires one ordinary data file and -one, with lines, scatter, or histogram"
    );
    ensure!(
        o.group_by.is_none() || o.missing != crate::input::Missing::Mean,
        "group-by does not impute across groups; choose error or gap"
    );
    ensure!(
        o.qc.lower.is_empty() == o.qc.upper.is_empty(),
        "supply both -yfiles_lo and -yfiles_hi"
    );
    ensure!(
        o.qc.pm.is_empty() || o.qc.lower.is_empty(),
        "choose symmetric widths or lower/upper bounds, not both"
    );
    ensure!(
        o.qc.lower.is_empty() || o.missing != crate::input::Missing::Mean,
        "supplied bounds cannot be mean-imputed; choose error or gap"
    );
    ensure!(
        !o.qc.use_pm || !o.qc.pm.is_empty() || !o.qc.lower.is_empty(),
        "ylim_use_pm requires supplied uncertainty widths or bounds"
    );
    ensure!(
        o.qc.box_view.is_empty() || o.qc.boxplot,
        "bplot_view requires boxplot_on"
    );
    ensure!(!o.qc.boxplot || !o.one, "boxplot_on needs separate panels");
    ensure!(
        !o.qc.margin_off || (o.caption.is_empty() && o.panel_captions.is_empty()),
        "margin_off cannot show captions"
    );
    ensure!(
        !o.qc.margin_off || (!o.qc.boxplot && o.panel_titles.is_empty()),
        "margin_off cannot show boxplots or panel titles"
    );
    ensure!(
        !o.diagonal || o.qc.ylimits.is_empty() || (o.qc.ylimits.len() == 1 && o.ylim.is_some()),
        "diagonal requires a complete single yaxis range or automatic limits"
    );
    ensure!(
        o.shades.is_empty() || (o.demo.is_none() && o.kind != "heat"),
        "-shade supports data-driven line, scatter, and histogram plots"
    );
    ensure!(
        o.shade_design.file.is_some() || o.shade_design == crate::shading::Design::default(),
        "shade-xmat settings require -shade-xmat FILE"
    );
    ensure!(
        o.shade_design.file.is_none() || (o.demo.is_none() && o.kind == "line"),
        "shade-xmat requires data-driven line plots"
    );
    let mut qc_drawing = o.qc.clone();
    qc_drawing.ylimits.clear();
    ensure!(
        qc_drawing == crate::qc::Options::default() || o.kind == "line",
        "QC annotations and scaling currently require line plots"
    );
    ensure!(
        o.demo.is_none()
            || (o.qc == crate::qc::Options::default() && o.missing == crate::input::Missing::Error),
        "demo cannot use QC or missing-data options"
    );
    if let Some(d) = &o.demo {
        ensure!(
            ["gallery", "line", "histogram", "scatter", "heat"].contains(&d.as_str()),
            "unknown demo '{d}'"
        );
        ensure!(
            !kind_set,
            "choose -demo NAME rather than a separate plot-type flag"
        );
        ensure!(
            !o.one
                && !o.sepscl
                && !o.reverse
                && o.xlabel.is_empty()
                && o.ylabels.is_empty()
                && o.legends.is_empty()
                && !o.legend_on
                && o.vlines.is_empty()
                && o.xlim.is_none()
                && o.ylim.is_none()
                && o.colors.is_empty()
                && o.dx == 1.
                && o.xzero == 0.,
            "demo scenes have fixed drawing settings; use examples/*.1D to customize plots"
        );
    }
    Ok(o)
}
