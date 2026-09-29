//! 2D pan/zoom camera mapping between world (network) coordinates and screen
//! (pixel) coordinates. Framework-agnostic: no egui/wgpu types appear here so
//! this can be reused by any renderer.

/// An axis-aligned bounding box in world space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bounds {
    pub min_x: f64,
    pub min_y: f64,
    pub max_x: f64,
    pub max_y: f64,
}

impl Bounds {
    pub fn from_points<I: IntoIterator<Item = (f64, f64)>>(points: I) -> Option<Self> {
        let mut it = points.into_iter();
        let (x0, y0) = it.next()?;
        let mut bounds = Bounds {
            min_x: x0,
            min_y: y0,
            max_x: x0,
            max_y: y0,
        };
        for (x, y) in it {
            bounds.min_x = bounds.min_x.min(x);
            bounds.min_y = bounds.min_y.min(y);
            bounds.max_x = bounds.max_x.max(x);
            bounds.max_y = bounds.max_y.max(y);
        }
        Some(bounds)
    }

    pub fn width(&self) -> f64 {
        self.max_x - self.min_x
    }
    pub fn height(&self) -> f64 {
        self.max_y - self.min_y
    }
    pub fn center(&self) -> (f64, f64) {
        (
            (self.min_x + self.max_x) / 2.0,
            (self.min_y + self.max_y) / 2.0,
        )
    }
    /// Expand the bounds by a fractional margin on each side (e.g. 0.1 = 10%).
    pub fn with_margin(&self, margin: f64) -> Self {
        let dx = (self.width() * margin).max(1.0);
        let dy = (self.height() * margin).max(1.0);
        Bounds {
            min_x: self.min_x - dx,
            min_y: self.min_y - dy,
            max_x: self.max_x + dx,
            max_y: self.max_y + dy,
        }
    }

    pub fn intersects(&self, other: &Bounds) -> bool {
        self.min_x <= other.max_x
            && self.max_x >= other.min_x
            && self.min_y <= other.max_y
            && self.max_y >= other.min_y
    }

    pub fn contains_point(&self, x: f64, y: f64) -> bool {
        x >= self.min_x && x <= self.max_x && y >= self.min_y && y <= self.max_y
    }
}

/// A simple affine 2D camera: `screen = (world - center) * zoom + screen_center`.
/// World Y increases upward (typical GIS/engineering convention); screen Y
/// increases downward, so the Y axis is flipped in the transform.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Camera {
    /// World-space point currently at the center of the viewport.
    pub center: (f64, f64),
    /// Screen pixels per world unit.
    pub zoom: f64,
    /// Size of the viewport in pixels (updated every frame by the caller).
    pub viewport: (f64, f64),
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            center: (0.0, 0.0),
            zoom: 1.0,
            viewport: (800.0, 600.0),
        }
    }
}

impl Camera {
    pub fn world_to_screen(&self, world: (f64, f64)) -> (f64, f64) {
        let (cx, cy) = self.center;
        let (vw, vh) = self.viewport;
        let x = (world.0 - cx) * self.zoom + vw / 2.0;
        let y = -(world.1 - cy) * self.zoom + vh / 2.0;
        (x, y)
    }

    pub fn screen_to_world(&self, screen: (f64, f64)) -> (f64, f64) {
        let (cx, cy) = self.center;
        let (vw, vh) = self.viewport;
        let x = (screen.0 - vw / 2.0) / self.zoom + cx;
        let y = -(screen.1 - vh / 2.0) / self.zoom + cy;
        (x, y)
    }

    /// Pan the camera by a screen-space pixel delta.
    pub fn pan_screen(&mut self, dx: f64, dy: f64) {
        self.center.0 -= dx / self.zoom;
        self.center.1 += dy / self.zoom;
    }

    /// Zoom by `factor` (>1 zooms in), keeping the world point currently under
    /// `screen_anchor` stationary on screen.
    pub fn zoom_at(&mut self, factor: f64, screen_anchor: (f64, f64)) {
        let world_anchor = self.screen_to_world(screen_anchor);
        self.zoom = (self.zoom * factor).clamp(1e-6, 1e6);
        let (vw, vh) = self.viewport;
        // Solve world_to_screen(world_anchor) == screen_anchor for `center`
        // directly under the new zoom, rather than zooming then panning:
        // panning by a screen delta divides by `zoom`, which after a zoom
        // change no longer corresponds to the pre-zoom screen delta.
        self.center.0 = world_anchor.0 - (screen_anchor.0 - vw / 2.0) / self.zoom;
        self.center.1 = world_anchor.1 + (screen_anchor.1 - vh / 2.0) / self.zoom;
    }

    /// The world-space bounds currently visible in the viewport.
    pub fn visible_bounds(&self) -> Bounds {
        let top_left = self.screen_to_world((0.0, 0.0));
        let bottom_right = self.screen_to_world(self.viewport);
        Bounds {
            min_x: top_left.0.min(bottom_right.0),
            min_y: top_left.1.min(bottom_right.1),
            max_x: top_left.0.max(bottom_right.0),
            max_y: top_left.1.max(bottom_right.1),
        }
    }

    /// Fit the camera to show `bounds` (with a small margin), keeping the
    /// current viewport size.
    pub fn fit_to_bounds(&mut self, bounds: Bounds) {
        let bounds = bounds.with_margin(0.08);
        self.center = bounds.center();
        let (vw, vh) = self.viewport;
        if bounds.width() <= 0.0 && bounds.height() <= 0.0 {
            self.zoom = 1.0;
            return;
        }
        let zoom_x = if bounds.width() > 0.0 {
            vw / bounds.width()
        } else {
            f64::MAX
        };
        let zoom_y = if bounds.height() > 0.0 {
            vh / bounds.height()
        } else {
            f64::MAX
        };
        self.zoom = zoom_x.min(zoom_y).clamp(1e-6, 1e6);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_world_screen() {
        let cam = Camera {
            center: (10.0, -5.0),
            zoom: 2.0,
            viewport: (400.0, 300.0),
        };
        let world = (12.5, 3.0);
        let screen = cam.world_to_screen(world);
        let back = cam.screen_to_world(screen);
        assert!((back.0 - world.0).abs() < 1e-9);
        assert!((back.1 - world.1).abs() < 1e-9);
    }

    #[test]
    fn zoom_at_keeps_anchor_fixed() {
        let mut cam = Camera::default();
        let anchor_screen = (100.0, 150.0);
        let anchor_world_before = cam.screen_to_world(anchor_screen);
        cam.zoom_at(2.0, anchor_screen);
        let anchor_world_after = cam.screen_to_world(anchor_screen);
        assert!((anchor_world_before.0 - anchor_world_after.0).abs() < 1e-9);
        assert!((anchor_world_before.1 - anchor_world_after.1).abs() < 1e-9);
    }

    #[test]
    fn fit_to_bounds_centers_camera() {
        let mut cam = Camera::default();
        let bounds = Bounds {
            min_x: 0.0,
            min_y: 0.0,
            max_x: 100.0,
            max_y: 50.0,
        };
        cam.fit_to_bounds(bounds);
        assert!((cam.center.0 - 50.0).abs() < 1e-6);
        assert!((cam.center.1 - 25.0).abs() < 1e-6);
    }
}
