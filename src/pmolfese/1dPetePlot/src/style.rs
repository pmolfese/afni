//! Presentation overrides. Source values, bounds, and scientific annotations stay untouched.
use crate::model::{Content, Plot};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};

pub const PALETTES: [&str; 3] = ["viridis", "Reds_and_Blues", "Reds_and_Blues_Inv"];
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct SeriesStyle {
    pub color: Option<String>,
    pub label: Option<String>,
    pub line_width: Option<f64>,
    pub marker_radius: Option<f64>,
}
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PlotStyle {
    pub series: Vec<SeriesStyle>,
    pub line_width: Option<f64>,
    pub marker_radius: Option<f64>,
    pub legend: Option<bool>,
    pub palette: Option<String>,
}
impl PlotStyle {
    pub fn validate(&self, plot: &Plot) -> Result<()> {
        ensure!(
            self.series.len() <= plot.style_series().len(),
            "style has more series than its graph"
        );
        let widths = matches!(plot.content, Content::Lines(_) | Content::Intervals { .. });
        let markers = widths || matches!(plot.content, Content::Scatter(_));
        let legends = matches!(
            plot.content,
            Content::Lines(_) | Content::Scatter(_) | Content::Histogram { .. }
        );
        ensure!(
            widths
                || (self.line_width.is_none()
                    && self.series.iter().all(|s| s.line_width.is_none())),
            "line width requires lines or intervals"
        );
        ensure!(
            markers
                || (self.marker_radius.is_none()
                    && self.series.iter().all(|s| s.marker_radius.is_none())),
            "marker size requires lines, scatter, or intervals"
        );
        ensure!(
            legends || (self.legend.is_none() && self.series.iter().all(|s| s.label.is_none())),
            "series legend requires lines, scatter, or histogram"
        );
        for width in self
            .line_width
            .iter()
            .chain(self.series.iter().filter_map(|s| s.line_width.as_ref()))
        {
            ensure!(
                width.is_finite() && (0.25..=12.).contains(width),
                "line width must be between 0.25 and 12"
            );
        }
        for radius in self
            .marker_radius
            .iter()
            .chain(self.series.iter().filter_map(|s| s.marker_radius.as_ref()))
        {
            ensure!(
                radius.is_finite() && (0.5..=20.).contains(radius),
                "marker radius must be between 0.5 and 20"
            );
        }
        for s in &self.series {
            if let Some(c) = &s.color {
                // Store canonical colors so renderer interpolation is always safe.
                ensure!(
                    crate::model::color(c).as_ref().is_ok_and(|v| v == c),
                    "style color must be canonical #rrggbb"
                );
            }
            if let Some(label) = &s.label {
                ensure!(label.len() <= 1024, "legend label exceeds 1024 bytes");
            }
        }
        if let Some(palette) = &self.palette {
            ensure!(
                PALETTES.contains(&palette.as_str()),
                "unknown style palette"
            );
            ensure!(
                matches!(plot.content, Content::Heat { .. }),
                "palette requires a heatmap"
            );
        }
        Ok(())
    }
    pub fn color<'a>(&'a self, index: usize, original: &'a str) -> &'a str {
        self.series
            .get(index)
            .and_then(|s| s.color.as_deref())
            .unwrap_or(original)
    }
    pub fn label<'a>(&'a self, index: usize, original: &'a str) -> &'a str {
        self.series
            .get(index)
            .and_then(|s| s.label.as_deref())
            .unwrap_or(original)
    }
    pub fn width(&self, index: usize, original: f64) -> f64 {
        self.series
            .get(index)
            .and_then(|s| s.line_width)
            .or(self.line_width)
            .unwrap_or(original)
    }
    pub fn radius(&self, index: usize, original: f64) -> f64 {
        self.series
            .get(index)
            .and_then(|s| s.marker_radius)
            .or(self.marker_radius)
            .unwrap_or(original)
    }
    pub fn cmap<'a>(&'a self, original: &'a str) -> &'a str {
        self.palette.as_deref().unwrap_or(original)
    }
}
impl Plot {
    /// Stable source series order, independent of user-facing legend names.
    pub fn style_series(&self) -> Vec<(&str, &str)> {
        match &self.content {
            Content::Lines(s) | Content::Scatter(s) | Content::Histogram { series: s, .. } => s
                .iter()
                .map(|s| (s.label.as_str(), s.color.as_str()))
                .collect(),
            Content::Intervals { color, .. } => vec![("Estimates and confidence intervals", color)],
            Content::Heat { .. } => vec![],
        }
    }
}
