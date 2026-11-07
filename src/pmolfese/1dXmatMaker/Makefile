# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

# Standalone wrapper: this does not add Rust to other AFNI programs' builds.
# Use PREFIX/DESTDIR explicitly when staging an install for validation.
CARGO ?= cargo
PREFIX ?= /usr/local
DESTDIR ?=
MAKER_TARGET_DIR ?= $(CURDIR)/target/make

.PHONY: all check install
all:
	$(CARGO) build --manifest-path "$(CURDIR)/Cargo.toml" --target-dir "$(MAKER_TARGET_DIR)" --release --frozen

check:
	$(CARGO) test --manifest-path "$(CURDIR)/Cargo.toml" --target-dir "$(MAKER_TARGET_DIR)" --frozen

install: all
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 "$(MAKER_TARGET_DIR)/release/1dXmatMaker" "$(DESTDIR)$(PREFIX)/bin/1dXmatMaker"

	install -d "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/licenses"
	install -m 644 README.md THIRD_PARTY.md "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/"
	# The waveform bank is embedded in the binary; install only its notices.
	install -d "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/data/glmsingle"
	install -m 644 data/glmsingle/LICENSE data/glmsingle/provenance.json "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/data/glmsingle/"
	install -d "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/examples"
	install -m 644 examples/*.tsv examples/*.1D "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/examples/"
	install -m 644 licenses/*.txt "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/licenses/"
	install -m 644 data/glmsingle/LICENSE "$(DESTDIR)$(PREFIX)/share/1dXmatMaker/licenses/GLMsingle.txt"
