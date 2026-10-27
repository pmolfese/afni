# Bundled plot font

`Ubuntu-Light.ttf` is the unmodified Ubuntu Light font distributed with
`epaint_default_fonts` 0.36.2 (egui). Its license is included as `UFL.txt`.
Font copyright: Canonical Ltd.; see the embedded font metadata and license.

The file exporter loads this bundled font only; it does not scan system fonts.
Exported SVG text is converted to vector paths, so recipients do not need the
font installed. This keeps appearance stable but makes labels outlines rather
than editable/searchable SVG text. The native UI uses egui's bundled fonts.
The plotting font does not cover every Unicode script; extended font selection
and missing-glyph diagnostics remain future work.

Upstream: https://github.com/emilk/egui/tree/0.36.2/crates/epaint_default_fonts/fonts
Font license: https://ubuntu.com/legal/font-licence

## Application icon

`peteplot.svg` is an original vector icon for this project: a white P over an
amber trace and teal plot axes. `peteplot.png` is its preview, and `peteplot.icns`
is used by the macOS bundle. No icon font, downloaded artwork, or runtime font
is needed. Regenerate on macOS with:

```sh
cargo run --locked --no-default-features --example render_icon
iconutil -c icns target/peteplot.iconset -o assets/peteplot.icns
```

The native window icon is rendered from the same embedded SVG.
