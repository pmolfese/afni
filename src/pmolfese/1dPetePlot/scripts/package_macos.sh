#!/bin/sh
# Build a self-contained local app. Does not install it into /Applications.
set -eu
PETEPLOT_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PETEPLOT_PROFILE=${1:-release}
case "$PETEPLOT_PROFILE" in
  release) cargo build --manifest-path "$PETEPLOT_ROOT/Cargo.toml" --locked --release ;;
  debug) cargo build --manifest-path "$PETEPLOT_ROOT/Cargo.toml" --locked ;;
  *) echo "usage: $0 [release|debug]" >&2; exit 1 ;;
esac
PETEPLOT_APP="$PETEPLOT_ROOT/target/1dPetePlot.app"
mkdir -p "$PETEPLOT_APP/Contents/MacOS"
mkdir -p "$PETEPLOT_APP/Contents/Resources"
cp "$PETEPLOT_ROOT/assets/peteplot.icns" "$PETEPLOT_APP/Contents/Resources/peteplot.icns"
# Replace the executable atomically so an already-running instance keeps its inode.
cp "$PETEPLOT_ROOT/target/$PETEPLOT_PROFILE/1dPetePlot" "$PETEPLOT_APP/Contents/MacOS/peteplot-bin.new"
mv -f "$PETEPLOT_APP/Contents/MacOS/peteplot-bin.new" "$PETEPLOT_APP/Contents/MacOS/peteplot-bin"
cat > "$PETEPLOT_APP/Contents/MacOS/1dPetePlot" <<'LAUNCHER'
#!/bin/sh
PETEPLOT_BIN_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
# Finder opens the empty import workspace. Pass -demo for the synthetic gallery.
exec "$PETEPLOT_BIN_DIR/peteplot-bin" "$@"
LAUNCHER
chmod +x "$PETEPLOT_APP/Contents/MacOS/1dPetePlot"
cat > "$PETEPLOT_APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
<key>CFBundleExecutable</key><string>1dPetePlot</string>
<key>CFBundleIdentifier</key><string>org.afni.1dPetePlot.prototype</string>
<key>CFBundleName</key><string>1dPetePlot</string>
<key>CFBundleDisplayName</key><string>1dPetePlot</string>
<key>CFBundleIconFile</key><string>peteplot</string>
<key>CFBundlePackageType</key><string>APPL</string>
<key>CFBundleVersion</key><string>0.1.0</string>
<key>NSHighResolutionCapable</key><true/>
</dict></plist>
PLIST
printf '%s\n' "$PETEPLOT_APP"
