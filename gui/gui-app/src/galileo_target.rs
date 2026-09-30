//! Bridges Galileo's `wgpu` renderer into an egui-drawn texture, sharing
//! the *same* `wgpu::Device`/`Queue` `eframe`'s own `egui-wgpu` backend
//! already created (see `app.rs::EditorApp::new`, which requires
//! `NativeOptions::renderer = Renderer::Wgpu` so `cc.wgpu_render_state` is
//! populated) - no second GPU context, no cross-API texture copies.
//!
//! The `Map` this targets currently has **no layers** - node/link
//! topology rendering moved to `instanced_renderer.rs`'s true GPU
//! instancing (one draw call per shape type, not a batched-bundle
//! CPU-tessellation approach; see its module doc for why). This module
//! and its `Map`/`WgpuRenderer` plumbing stay in place regardless,
//! because a `.pmtiles`/vector-tile basemap layer is still a planned,
//! Galileo-appropriate use (real GIS features, MapLibre-style rendering -
//! exactly what Galileo *is* built for, unlike the fixed small shape
//! catalog `instanced_renderer.rs` now owns): add that layer to the `Map`
//! built in `app.rs::EditorApp::new` when it exists, and it composites as
//! the backdrop underneath the instanced topology `canvas.rs` paints on
//! top of this module's output - no other changes needed here.
//!
//! Deliberately does **not** use `galileo_egui::EguiMap`/`EguiMapState`.
//! That widget is a real, supported integration path, but it owns pointer
//! input itself (`EventProcessor` + `MapController`, feeding
//! `UserEventHandler`s) via its own `ui.allocate_exact_size` call - which
//! would fight `canvas.rs`'s existing `ui.allocate_rect` + gesture handling
//! for the same screen region (both trying to own the same drag/click).
//! Using `galileo::render::WgpuRenderer` directly instead - the same calls
//! `EguiMapState::render()` makes internally - keeps drawing and input
//! fully decoupled: `canvas.rs` still allocates the rect and runs its own
//! Select/AddJunction/AddLink/box-select handling completely unchanged;
//! this module only turns the current `Map` into a texture and hands back
//! an `egui::TextureId` for `canvas.rs` to paint.

use eframe::egui;
use galileo::Map;
use galileo::render::WgpuRenderer;
use galileo_types::cartesian::Size;

/// Owns the Galileo map + its offscreen render target, and the egui
/// texture that target is registered as.
///
/// `renderer`/`size`/`last_rendered_key`/`paint` (below) are currently
/// unused - `canvas.rs` no longer calls `paint` (see its basemap-backdrop
/// comment) since painting an empty map actively overwrote the real
/// canvas background. Kept, not removed, for when a real basemap layer
/// makes calling it again worthwhile.
#[allow(dead_code)]
pub struct GalileoTarget {
    renderer: WgpuRenderer,
    map: Map,
    texture_id: Option<egui::TextureId>,
    size: (u32, u32),
    /// `(map center x/y, resolution, content revision)` from the last call
    /// that actually rendered - see `paint`'s doc comment. `f64`s are
    /// compared for bit-exact equality deliberately: they're read straight
    /// from `Camera` each frame, so an unchanged camera reliably produces
    /// the identical value, not a value merely close to it.
    last_rendered_key: Option<(f64, f64, f64, u64)>,
}

impl GalileoTarget {
    /// `render_state` is `cc.wgpu_render_state.clone().expect(...)` from
    /// `eframe::CreationContext` - see `app.rs`.
    pub fn new(render_state: &egui_wgpu::RenderState, map: Map) -> Self {
        let device = render_state.device.clone();
        let queue = render_state.queue.clone();
        // 1x1 real target; `resize` (called every frame from `paint`)
        // grows it to the canvas rect's actual pixel size before the first
        // real render, so this initial size is never actually drawn to.
        let renderer = WgpuRenderer::new_with_device_and_texture(device, queue, Size::new(1, 1));

        Self {
            renderer,
            map,
            texture_id: None,
            size: (1, 1),
            last_rendered_key: None,
        }
    }

    pub fn map(&mut self) -> &mut Map {
        &mut self.map
    }

    /// Renders the current `Map` state at `size_px` and paints it into
    /// `rect` via `painter` - unless neither the view (`map_center_x/y`,
    /// `resolution`) nor the content (`content_revision`, from
    /// `NetworkLayers::revision`) has changed since the last call that
    /// actually rendered, in which case the previous frame's texture is
    /// re-blitted as-is and no GPU work happens.
    ///
    /// That check is load-bearing, not an optimization on top of an
    /// already-fine baseline: `canvas.rs` calls this once per egui frame,
    /// and egui repaints continuously while the pointer moves (hover
    /// highlighting, tooltips, ...) - including when the pointer is over
    /// a totally unrelated panel/button, nowhere near the canvas. Without
    /// this, every one of those repaints would re-render the entire
    /// network (rasterizing however many hundred thousand already-
    /// tessellated primitives are visible) for no reason, which is a much
    /// larger, unconditional per-frame cost than anything
    /// `NetworkLayers::sync`'s dirty-tracking governs - see its own doc
    /// comment for why *tessellation* itself is properly cached by Galileo
    /// and isn't the issue here; re-rendering/re-presenting an unchanged
    /// scene is a separate cost this module has to avoid on its own.
    #[allow(dead_code)]
    pub fn paint(
        &mut self,
        render_state: &egui_wgpu::RenderState,
        painter: &egui::Painter,
        rect: egui::Rect,
        size_px: (u32, u32),
        map_center: (f64, f64),
        resolution: f64,
        content_revision: u64,
    ) {
        if size_px.0 == 0 || size_px.1 == 0 {
            return;
        }

        let resized = size_px != self.size;
        if resized {
            self.renderer.resize(Size::new(size_px.0, size_px.1));
            self.size = size_px;
            // The old egui texture id (if any) pointed at a `TextureView`
            // sized for the *previous* resolution; `wgpu` doesn't let a
            // view outlive the texture it was created from being replaced
            // by `resize`, so re-register under a fresh id rather than
            // try to update the stale one in place.
            if let Some(id) = self.texture_id.take() {
                render_state.renderer.write().free_texture(&id);
            }
        }

        let key = (map_center.0, map_center.1, resolution, content_revision);
        let needs_render = resized || self.texture_id.is_none() || self.last_rendered_key != Some(key);

        if needs_render {
            self.map.set_size(Size::new(size_px.0 as f64, size_px.1 as f64));

            let Some(view) = self.renderer.get_target_texture_view() else {
                return;
            };
            self.renderer.render_to_texture_view(&self.map, &view);
            self.last_rendered_key = Some(key);

            if self.texture_id.is_none() {
                let id = render_state.renderer.write().register_native_texture(
                    &render_state.device,
                    &view,
                    wgpu::FilterMode::Linear,
                );
                self.texture_id = Some(id);
            }
        }

        let Some(texture_id) = self.texture_id else {
            return;
        };

        painter.image(
            texture_id,
            rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            egui::Color32::WHITE,
        );
    }
}

impl Drop for GalileoTarget {
    fn drop(&mut self) {
        // Best-effort: `render_state` isn't available in `Drop`, so the
        // texture is simply abandoned; `egui-wgpu` doesn't require
        // `free_texture` to be called for correctness, only to reclaim the
        // id promptly - fine for an editor with one long-lived canvas.
        let _ = &self.texture_id;
    }
}
