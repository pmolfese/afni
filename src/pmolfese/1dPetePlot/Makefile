# Standalone prototype wrapper. AFNI's top-level build is not modified yet.
CARGO ?= cargo
PREFIX ?= /usr/local
GUI ?= 1
CARGO_FLAGS := --release --locked
ifeq ($(GUI),0)
CARGO_FLAGS += --no-default-features
endif

.PHONY: all test gallery app install clean
all:
	$(CARGO) build $(CARGO_FLAGS)

test:
	$(CARGO) test --locked
	$(CARGO) test --locked --no-default-features

gallery: all
	./scripts/render_gallery.sh

app:
	./scripts/package_macos.sh release

install: all
	install -d "$(DESTDIR)$(PREFIX)/bin"
	install -m 755 target/release/1dPetePlot "$(DESTDIR)$(PREFIX)/bin/1dPetePlot"

clean:
	$(CARGO) clean
