//! macOS owns a global AppKit menu. Linux uses the same actions in-window.
use eframe::egui;
use muda::{
    Menu, MenuItem, PredefinedMenuItem as Native, Submenu,
    accelerator::{Accelerator, Code, Modifiers},
};

pub struct NativeMenu {
    _menu: Menu,
    fullscreen: MenuItem,
    events: std::sync::mpsc::Receiver<muda::MenuEvent>,
}
impl NativeMenu {
    pub fn new(ctx: &egui::Context) -> anyhow::Result<Self> {
        let menu = Menu::new();
        let app = Submenu::with_items(
            "1dPetePlot",
            true,
            &[
                &Native::about(
                    Some("About 1dPetePlot"),
                    Some(muda::AboutMetadata {
                        name: Some("1dPetePlot".into()),
                        version: Some(env!("CARGO_PKG_VERSION").into()),
                        comments: Some("Native plotting for AFNI".into()),
                        ..Default::default()
                    }),
                ),
                &Native::separator(),
                &Native::hide(None),
                &Native::hide_others(None),
                &Native::show_all(None),
                &Native::separator(),
                &MenuItem::with_id(
                    "quit",
                    "Quit 1dPetePlot",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::KeyQ)),
                ),
            ],
        )?;
        let save = MenuItem::with_id(
            "save",
            "Save Image…",
            true,
            Some(Accelerator::new(Modifiers::META, Code::KeyS)),
        );
        let file = Submenu::with_items(
            "File",
            true,
            &[
                &MenuItem::with_id(
                    "data-open",
                    "Open Data…",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::KeyI)),
                ),
                &MenuItem::with_id(
                    "project-open",
                    "Open Project…",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::KeyO)),
                ),
                &MenuItem::with_id(
                    "project-save",
                    "Save Project…",
                    true,
                    Some(Accelerator::new(
                        Modifiers::META | Modifiers::SHIFT,
                        Code::KeyS,
                    )),
                ),
                &MenuItem::with_id("project-save-as", "Save Project As…", true, None),
                &Native::separator(),
                &save,
                &Native::separator(),
                &MenuItem::with_id(
                    "quit",
                    "Close Window",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::KeyW)),
                ),
            ],
        )?;
        let edit = Submenu::with_items(
            "Edit",
            true,
            &[
                &MenuItem::with_id(
                    "undo",
                    "Undo",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::KeyZ)),
                ),
                &MenuItem::with_id(
                    "redo",
                    "Redo",
                    true,
                    Some(Accelerator::new(
                        Modifiers::META | Modifiers::SHIFT,
                        Code::KeyZ,
                    )),
                ),
                &Native::separator(),
                &Native::cut(None),
                &Native::copy(None),
                &Native::paste(None),
                &Native::select_all(None),
                &Native::separator(),
                &MenuItem::with_id("copy", "Copy Selected Graph", true, None),
                &MenuItem::with_id("copy-all", "Copy All Graphs", true, None),
            ],
        )?;
        let fullscreen = MenuItem::with_id(
            "fullscreen",
            "Enter Full Screen",
            true,
            Some(Accelerator::new(
                Modifiers::META | Modifiers::CONTROL,
                Code::KeyF,
            )),
        );
        let view = Submenu::with_items(
            "View",
            true,
            &[
                &MenuItem::with_id(
                    "reset",
                    "Reset View",
                    true,
                    Some(Accelerator::new(Modifiers::META, Code::Digit0)),
                ),
                &MenuItem::with_id("overview", "Show All Graphs", true, None),
                &Native::separator(),
                &fullscreen,
            ],
        )?;
        let window = Submenu::with_items(
            "Window",
            true,
            &[&Native::minimize(None), &Native::maximize(None)],
        )?;
        menu.append_items(&[&app, &file, &edit, &view, &window])?;
        menu.init_for_nsapp();
        window.set_as_windows_menu_for_nsapp();
        let (tx, events) = std::sync::mpsc::channel();
        let ctx = ctx.clone();
        muda::MenuEvent::set_event_handler(Some(move |event| {
            let _ = tx.send(event);
            ctx.request_repaint();
        }));
        Ok(Self {
            _menu: menu,
            fullscreen,
            events,
        })
    }
    pub fn next(&self) -> Option<String> {
        self.events.try_recv().ok().map(|e| e.id.0)
    }
    pub fn set_fullscreen(&self, fullscreen: bool) {
        let title = if fullscreen {
            "Exit Full Screen"
        } else {
            "Enter Full Screen"
        };
        if self.fullscreen.text() != title {
            self.fullscreen.set_text(title);
        }
    }
}
