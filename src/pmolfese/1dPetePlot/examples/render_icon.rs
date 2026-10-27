// Regenerate the iconset from the editable vector source; no desktop required.
use resvg::{tiny_skia, usvg};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let dir = root.join("target/peteplot.iconset");
    std::fs::create_dir_all(&dir)?;
    let tree = usvg::Tree::from_str(
        include_str!("../assets/peteplot.svg"),
        &usvg::Options::default(),
    )?;
    for size in [16, 32, 128, 256, 512] {
        for scale in [1, 2] {
            let n = size * scale;
            let mut pix = tiny_skia::Pixmap::new(n, n).unwrap();
            resvg::render(
                &tree,
                tiny_skia::Transform::from_scale(n as f32 / 256., n as f32 / 256.),
                &mut pix.as_mut(),
            );
            let suffix = if scale == 2 { "@2x" } else { "" };
            pix.save_png(dir.join(format!("icon_{size}x{size}{suffix}.png")))?;
            if n == 256 && scale == 1 {
                pix.save_png(root.join("assets/peteplot.png"))?;
            }
        }
    }
    Ok(())
}
