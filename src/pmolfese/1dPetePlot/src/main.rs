#[cfg(feature = "gui")]
mod attachments;
mod cli;
mod document;
mod grouping;
#[cfg(feature = "gui")]
mod gui;
#[cfg(feature = "gui")]
mod import;
#[cfg(feature = "gui")]
mod import_ui;
mod input;
#[cfg(all(feature = "gui", target_os = "macos"))]
mod menus;
mod model;
mod qc;
mod rdm;
mod render;
mod rsa;
mod rsa_extra;
mod rsa_long;
mod shading;
mod style;

use anyhow::{Context, Result};
fn main() {
    if let Err(e) = run() {
        eprintln!("** 1dPetePlot: {e:#}");
        std::process::exit(1);
    }
}
fn run() -> Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args
        .iter()
        .any(|a| ["-help", "--help", "-h"].contains(&a.as_str()))
    {
        print!("{}", cli::HELP);
        return Ok(());
    }
    if args == ["-version"] {
        println!(
            "1dPetePlot {} (prototype; gui={})",
            env!("CARGO_PKG_VERSION"),
            cfg!(feature = "gui")
        );
        return Ok(());
    }
    if args.is_empty() {
        #[cfg(feature = "gui")]
        return gui::run(
            model::Figure::empty(),
            cli::Options::default(),
            render::Renderer::new(),
            None,
            None,
            Some(Ok(document::Origin::empty())),
        );
        #[cfg(not(feature = "gui"))]
        anyhow::bail!("this build has no native-window support; supply input and -prefix plot.png");
    }
    let loaded = if args.first().is_some_and(|a| a == "-recipe") {
        let path = args.get(1).context("-recipe requires a project filename")?;
        Some(document::load(std::path::Path::new(path))?)
    } else {
        None
    };
    let mut o = if let Some(project) = &loaded {
        let overrides = &args[2..];
        let mut i = 0;
        while i < overrides.len() {
            match overrides[i].as_str() {
                "-show" | "-overwrite" => i += 1,
                "-prefix" | "-png" | "-svg" | "-width" | "-height" | "-dpi" => i += 2,
                "-figsize" => i += 3,
                flag => anyhow::bail!(
                    "{flag} cannot override a project; edit the project or use output/size options"
                ),
            }
        }
        let mut base = vec!["1D: 1 2".to_string()];
        if project
            .figure
            .plots
            .first()
            .is_some_and(|p| p.qc.margin_off)
        {
            base.push("-margin_off".into());
        }
        let mut options = cli::parse(&[base, overrides.to_vec()].concat())?;
        if options.auto_size {
            options.width = project.state.width;
            options.height = project.state.height;
        } else if !overrides.iter().any(|a| a == "-figsize" || a == "-dpi") {
            if !overrides.iter().any(|a| a == "-width") {
                options.width = project.state.width;
            }
            if !overrides.iter().any(|a| a == "-height") {
                options.height = project.state.height;
            }
        }
        options.auto_size = false;
        anyhow::ensure!(
            options.width as u64 * options.height as u64 <= 32_000_000,
            "figure exceeds 32 million pixels"
        );
        options
    } else {
        cli::parse(&args)?
    };
    // Fingerprint before reading plot data; later saves verify the same bytes.
    #[cfg(feature = "gui")]
    let origin = if loaded.is_none() && (o.output.is_none() || o.show) {
        Some(document::Origin::capture(&o))
    } else {
        None
    };
    let fig = if let Some(project) = &loaded {
        project.figure.clone()
    } else {
        model::from_options(&o).context("preparing plot")?
    };
    if o.auto_size && !o.qc.margin_off {
        let rows = fig.plots.len().div_ceil(fig.columns);
        o.height = o
            .height
            .max(rsa::minimum_height(&fig, o.width))
            .max(((165. + rows as f64 * 270.) * o.width as f64 / 1200.).ceil() as u32);
        anyhow::ensure!(
            o.height <= 8192 && o.width as u64 * o.height as u64 <= 32_000_000,
            "automatic layout exceeds image limits; use -columns, smaller -width, or fewer panels"
        );
    }
    let r = render::Renderer::new();
    let views = if let Some(project) = &loaded {
        if project.state.full_plot {
            fig.plots.iter().map(|p| p.bounds).collect::<Vec<_>>()
        } else {
            project.state.views.clone()
        }
    } else {
        fig.plots.iter().map(|p| p.bounds).collect::<Vec<_>>()
    };
    if let Some(path) = &o.output {
        r.export(
            &fig,
            &views,
            o.width,
            o.height,
            std::path::Path::new(path),
            o.overwrite,
        )?;
        println!("++ Saved {} × {}: {path}", o.width, o.height);
    }
    if o.output.is_none() || o.show {
        #[cfg(feature = "gui")]
        {
            gui::run(
                fig,
                o,
                r,
                loaded,
                args.get(1)
                    .filter(|_| args.first().is_some_and(|a| a == "-recipe"))
                    .map(std::path::PathBuf::from),
                origin,
            )?;
        }
        #[cfg(not(feature = "gui"))]
        anyhow::bail!(
            "this build has no native-window support; use -prefix plot.png or rebuild with GUI features"
        );
    }
    Ok(())
}
