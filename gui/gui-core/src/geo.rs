//! Web Mercator projection math and the network -> basemap georeferencing
//! transform.
//!
//! Pure math, no `egui`/`wgpu`/tile-format dependency (kept framework- and
//! basemap-library-agnostic, like the rest of this crate); `gui-app`'s
//! basemap layer uses this to place decoded MVT tile geometry (which is in
//! Web Mercator meters, the EPSG:3857 coordinate system essentially every
//! vector basemap - including Protomaps - uses) into the network's own
//! local/INP coordinate space, where the existing `Camera` already knows
//! how to draw it.
//!
//! `[COORDINATES]` in an INP file are not necessarily lon/lat (many real
//! networks use an arbitrary local or national grid - see the plan's
//! assumption 0), so this can't be a fixed formula: `GeoTransform::fit` is
//! the "pick 2+ known points, compute an affine transform" calibration
//! step the plan calls out as unavoidable.

use crate::camera::Bounds;

/// Half the circumference of the Web Mercator projection's square, in
/// meters (`EARTH_RADIUS * PI`). Web Mercator X and Y both range over
/// `[-ORIGIN_SHIFT, ORIGIN_SHIFT]`.
pub const ORIGIN_SHIFT: f64 = 20_037_508.342_789_244;
const EARTH_RADIUS: f64 = 6_378_137.0;
/// Web Mercator is undefined at the poles; real-world basemaps (including
/// Protomaps) clip latitude to this bound, matching the standard
/// projection's usable range.
const MAX_LATITUDE: f64 = 85.051_128_78;

/// Projects WGS84 longitude/latitude (degrees) to Web Mercator meters
/// (EPSG:3857) - the coordinate system tile `(z, x, y)` bounds are defined
/// in.
#[must_use]
pub fn lonlat_to_mercator(lon: f64, lat: f64) -> (f64, f64) {
    let x = lon.to_radians() * EARTH_RADIUS;
    let lat = lat.clamp(-MAX_LATITUDE, MAX_LATITUDE);
    let y = EARTH_RADIUS * (std::f64::consts::FRAC_PI_4 + lat.to_radians() / 2.0).tan().ln();
    (x, y)
}

/// The Web Mercator meter bounds of slippy-map tile `(z, x, y)` (XYZ tiling
/// scheme: `y` increases southward from the grid's top-left, matching both
/// `pmtiles::TileCoord` and the Protomaps/MVT convention).
#[must_use]
pub fn tile_bounds_mercator(z: u8, x: u32, y: u32) -> Bounds {
    let n = (1u64 << z) as f64;
    let tile_size = 2.0 * ORIGIN_SHIFT / n;
    let min_x = -ORIGIN_SHIFT + x as f64 * tile_size;
    let max_y = ORIGIN_SHIFT - y as f64 * tile_size;
    Bounds {
        min_x,
        min_y: max_y - tile_size,
        max_x: min_x + tile_size,
        max_y,
    }
}

/// A best-fit 2D similarity transform (uniform scale + rotation +
/// translation; no shear or reflection) mapping the network's own local
/// world/INP coordinates to Web Mercator meters.
///
/// Fit from 2+ `(network point, known lon/lat)` correspondences via the
/// standard closed-form least-squares similarity (2D orthogonal Procrustes
/// without reflection) - exact for exactly 2 points, a least-squares best
/// fit for more. This is deliberately a *global* linear approximation, not
/// a true reprojection: adequate for a single network's extent (city/
/// utility-district scale, what this editor targets), where the true
/// (non-linear) Mercator projection is locally close to affine anyway.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoTransform {
    scale: f64,
    cos_theta: f64,
    sin_theta: f64,
    tx: f64,
    ty: f64,
}

impl GeoTransform {
    /// Fits a transform from `control_points`, each a `(network_xy, (lon,
    /// lat))` pair. Returns `None` if there are fewer than 2 points, or the
    /// network points are degenerate (coincident, leaving no well-defined
    /// scale/rotation) - the caller should surface that as "add another,
    /// distinct point" rather than silently keeping a stale transform.
    #[must_use]
    pub fn fit(control_points: &[((f64, f64), (f64, f64))]) -> Option<Self> {
        if control_points.len() < 2 {
            return None;
        }
        let points: Vec<((f64, f64), (f64, f64))> = control_points
            .iter()
            .map(|&(world, (lon, lat))| (world, lonlat_to_mercator(lon, lat)))
            .collect();

        let n = points.len() as f64;
        let (wcx, wcy) = points
            .iter()
            .fold((0.0, 0.0), |(ax, ay), &((x, y), _)| (ax + x, ay + y));
        let (wcx, wcy) = (wcx / n, wcy / n);
        let (mcx, mcy) = points
            .iter()
            .fold((0.0, 0.0), |(ax, ay), &(_, (x, y))| (ax + x, ay + y));
        let (mcx, mcy) = (mcx / n, mcy / n);

        let mut num_cos = 0.0;
        let mut num_sin = 0.0;
        let mut denom = 0.0;
        for &((wx, wy), (mx, my)) in &points {
            let (dwx, dwy) = (wx - wcx, wy - wcy);
            let (dmx, dmy) = (mx - mcx, my - mcy);
            num_cos += dwx * dmx + dwy * dmy;
            num_sin += dwx * dmy - dwy * dmx;
            denom += dwx * dwx + dwy * dwy;
        }
        if denom < 1e-9 {
            return None;
        }
        let scale_cos = num_cos / denom;
        let scale_sin = num_sin / denom;
        let scale = scale_cos.hypot(scale_sin);
        if scale < 1e-12 {
            return None;
        }
        let cos_theta = scale_cos / scale;
        let sin_theta = scale_sin / scale;
        let tx = mcx - scale * (cos_theta * wcx - sin_theta * wcy);
        let ty = mcy - scale * (sin_theta * wcx + cos_theta * wcy);
        Some(Self { scale, cos_theta, sin_theta, tx, ty })
    }

    /// Maps a network/world point to Web Mercator meters.
    #[must_use]
    pub fn to_mercator(&self, world: (f64, f64)) -> (f64, f64) {
        let (x, y) = world;
        (
            self.scale * (self.cos_theta * x - self.sin_theta * y) + self.tx,
            self.scale * (self.sin_theta * x + self.cos_theta * y) + self.ty,
        )
    }

    /// Maps a Web Mercator point back to network/world coordinates - the
    /// inverse of `to_mercator`, used to place decoded tile geometry (in
    /// mercator meters) into world space for drawing with the existing
    /// `Camera`.
    #[must_use]
    pub fn to_world(&self, mercator: (f64, f64)) -> (f64, f64) {
        let (dx, dy) = (mercator.0 - self.tx, mercator.1 - self.ty);
        (
            (self.cos_theta * dx + self.sin_theta * dy) / self.scale,
            (-self.sin_theta * dx + self.cos_theta * dy) / self.scale,
        )
    }

    /// Mercator meters per world/network unit - used to pick which tile
    /// zoom level's on-disk resolution best matches the current camera
    /// zoom.
    #[must_use]
    pub fn scale(&self) -> f64 {
        self.scale
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Inverse of `lonlat_to_mercator`, only precise enough for this
    /// test's round-trip check (not part of the public API - real code
    /// never needs mercator -> lon/lat).
    fn mercator_to_lonlat_for_test(mercator: (f64, f64)) -> (f64, f64) {
        let (mx, my) = mercator;
        let lon = (mx / EARTH_RADIUS).to_degrees();
        let lat =
            (2.0 * (my / EARTH_RADIUS).exp().atan() - std::f64::consts::FRAC_PI_2).to_degrees();
        (lon, lat)
    }

    #[test]
    fn fits_exact_transform_from_two_points() {
        let world_a = (0.0, 0.0);
        let world_b = (10.0, 0.0);
        let mercator_a = (1000.0, 2000.0);
        // world +x maps to mercator +y: a 90-degree rotation.
        let mercator_b = (1000.0, 2010.0);

        let lonlat_a = mercator_to_lonlat_for_test(mercator_a);
        let lonlat_b = mercator_to_lonlat_for_test(mercator_b);

        let transform = GeoTransform::fit(&[(world_a, lonlat_a), (world_b, lonlat_b)]).unwrap();

        let got_a = transform.to_mercator(world_a);
        let got_b = transform.to_mercator(world_b);
        assert!((got_a.0 - mercator_a.0).abs() < 1e-6);
        assert!((got_a.1 - mercator_a.1).abs() < 1e-6);
        assert!((got_b.0 - mercator_b.0).abs() < 1e-6);
        assert!((got_b.1 - mercator_b.1).abs() < 1e-6);
    }

    #[test]
    fn to_world_inverts_to_mercator() {
        let transform = GeoTransform {
            scale: 2.5,
            cos_theta: 0.6,
            sin_theta: 0.8, // cos^2 + sin^2 == 1
            tx: 500.0,
            ty: -300.0,
        };
        let world = (12.5, -7.0);
        let mercator = transform.to_mercator(world);
        let back = transform.to_world(mercator);
        assert!((back.0 - world.0).abs() < 1e-9);
        assert!((back.1 - world.1).abs() < 1e-9);
    }

    #[test]
    fn fit_rejects_fewer_than_two_points() {
        assert!(GeoTransform::fit(&[]).is_none());
        assert!(GeoTransform::fit(&[((0.0, 0.0), (0.0, 0.0))]).is_none());
    }

    #[test]
    fn fit_rejects_coincident_points() {
        let p = ((5.0, 5.0), (4.0, 52.0));
        assert!(GeoTransform::fit(&[p, p]).is_none());
    }

    #[test]
    fn tile_bounds_cover_the_full_mercator_square_at_zoom_zero() {
        let bounds = tile_bounds_mercator(0, 0, 0);
        assert!((bounds.min_x + ORIGIN_SHIFT).abs() < 1e-6);
        assert!((bounds.max_x - ORIGIN_SHIFT).abs() < 1e-6);
        assert!((bounds.min_y + ORIGIN_SHIFT).abs() < 1e-6);
        assert!((bounds.max_y - ORIGIN_SHIFT).abs() < 1e-6);
    }

    #[test]
    fn tile_bounds_tile_at_z1_00_is_northwest_quadrant() {
        // At zoom 1, tile (0,0) is the northwest quarter of the world.
        let bounds = tile_bounds_mercator(1, 0, 0);
        assert!((bounds.min_x + ORIGIN_SHIFT).abs() < 1e-6);
        assert!((bounds.max_x).abs() < 1e-6);
        assert!((bounds.max_y - ORIGIN_SHIFT).abs() < 1e-6);
        assert!((bounds.min_y).abs() < 1e-6);
    }
}
