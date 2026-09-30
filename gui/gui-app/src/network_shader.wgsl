// Instanced network renderer: two draw calls total for the whole network
// (regardless of node/link count) - one instanced draw for "billboards"
// (node symbols + direction-arrow triangles, both fixed-pixel-size shapes
// oriented per-instance) and one for links (thick lines, screen-space
// width). See `instanced_renderer.rs`'s module doc for why this replaced
// the earlier Galileo-`FeatureLayer`-based approach: this is true GPU
// instancing (one small per-shape vertex buffer, reused via per-instance
// data, one draw call per shape *type* rather than per few-thousand
// features) - "render the whole network in one pass" in the way a game
// engine would, not a batched-bundle CPU-tessellation approach.
//
// Both shaders share this one module (and the camera uniform) so there's
// only one place the world -> screen -> NDC math lives.
//
// Shapes are drawn via a signed-distance-ish test in the fragment shader
// (a quad per instance, discarded/shaped down to a circle/triangle/thick
// line) rather than exactly-fitted geometry, so edges are aliased by
// default. `fwidth()` (a screen-space derivative, giving "how much does
// this value change per pixel") is used throughout to soften those edges
// over roughly a one-pixel band, at any zoom level - the standard
// SDF-antialiasing technique, and the only kind of antialiasing available
// here: this draws directly into egui's own render pass via
// `egui_wgpu::CallbackTrait::paint` (see `instanced_renderer.rs`), whose
// sample count is fixed by egui itself (1x, no MSAA) - a pipeline used
// within that pass must match it, so true MSAA isn't an option to just
// turn on for these two pipelines specifically.

struct Camera {
    // World-space point at the center of the viewport - matches
    // `gui_core::camera::Camera::world_to_screen`'s exact convention:
    // world Y increases *up*, screen Y increases *down* (flipped below),
    // and `center` really is the screen-center world point, not a corner.
    // Both this and every instance's own world-space fields are already
    // relative to `InstancedRenderer`'s floating origin (see its module
    // doc) - not raw INP/model coordinates - so this subtraction and the
    // rest of the pipeline stay in small, well-conditioned `f32` numbers
    // even for real-world-scale networks (e.g. Dutch RD coordinates, in
    // the hundreds of thousands of meters).
    center: vec2<f32>,
    zoom: f32,       // screen pixels per world unit
    _pad0: f32,
    viewport: vec2<f32>, // pixels
    _pad1: vec2<f32>,
};

@group(0) @binding(0) var<uniform> camera: Camera;

fn world_to_screen(world: vec2<f32>) -> vec2<f32> {
    let rel = world - camera.center;
    return vec2<f32>(rel.x, -rel.y) * camera.zoom + camera.viewport * 0.5;
}

fn screen_to_ndc(screen: vec2<f32>) -> vec4<f32> {
    let ndc = vec2<f32>(
        screen.x / camera.viewport.x * 2.0 - 1.0,
        1.0 - screen.y / camera.viewport.y * 2.0,
    );
    return vec4<f32>(ndc, 0.0, 1.0);
}

/// Soft edge over roughly a one-screen-pixel band around `dist == 0`
/// (negative = inside, positive = outside) - the shared antialiasing
/// primitive every shape test below uses.
fn aa_alpha(dist: f32) -> f32 {
    let width = max(fwidth(dist), 1e-6);
    return 1.0 - smoothstep(-width, width, dist);
}

// --- billboards: node symbols + direction arrows ---
//
// A fixed-pixel-size quad per instance, oriented by `direction` (already
// in *world* Y-up convention; flipped to screen convention below before
// use) and shaped by `kind` in the fragment shader via a cheap
// signed-distance-ish test against the *unrotated* local quad coordinate
// (the rotation already happened geometrically in the vertex shader, so
// the fragment shader can always test against one fixed-orientation
// shape regardless of `direction`).

struct BillboardVertexIn {
    @location(0) local_pos: vec2<f32>, // unit quad corner, [-1, 1]
};

struct BillboardInstanceIn {
    @location(1) center: vec2<f32>,     // world space (origin-relative - see `Camera`)
    @location(2) direction: vec2<f32>,  // world-space unit vector; (0,1) if orientation doesn't matter (circle/square)
    @location(3) size: f32,             // pixel diameter
    @location(4) color: vec4<f32>,
    @location(5) kind: u32,             // 0 = circle, 1 = square, 2 = triangle (reservoir), 3 = triangle (direction arrow)
    @location(6) world_length: f32,     // kind == 3 only: the link's route length, for the screen-space visibility cull below
};

struct BillboardVertexOut {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) local_pos: vec2<f32>,
    @location(1) color: vec4<f32>,
    @location(2) kind: u32,
};

/// Below this on-screen pipe length, a direction arrow is more clutter
/// than signal (it visually swallows the junction billboards at either
/// end - the reported "arrows overlap junctions when zoomed out") -
/// hidden (zero size) rather than drawn undersized.
const MIN_ARROW_SCREEN_LENGTH_PX: f32 = 28.0;

@vertex
fn vs_billboard(vin: BillboardVertexIn, iin: BillboardInstanceIn) -> BillboardVertexOut {
    // World Y-up tangent -> screen Y-down convention, matching
    // `world_to_screen`'s own flip, so a triangle "pointing" along a
    // link's world tangent actually points that way on screen.
    //
    // The fragment shader's shape test (`fs_billboard`) is defined with
    // the triangle's tip at local (0, 1) - i.e. its "forward" is local
    // +Y, not +X - so this needs the rotation that maps (0, 1) onto
    // `dir_screen` specifically (not the more common "rotate +X onto
    // dir" complex-multiplication form, which would point 90 degrees off
    // for every non-default direction).
    let dir_screen = normalize(vec2<f32>(iin.direction.x, -iin.direction.y));
    let rotated = vec2<f32>(
        vin.local_pos.x * dir_screen.y + vin.local_pos.y * dir_screen.x,
        -vin.local_pos.x * dir_screen.x + vin.local_pos.y * dir_screen.y,
    );

    var size = iin.size;
    if (iin.kind == 3u) {
        let screen_length = iin.world_length * camera.zoom;
        if (screen_length < MIN_ARROW_SCREEN_LENGTH_PX) {
            size = 0.0;
        }
    }

    let screen_center = world_to_screen(iin.center);
    let screen_pos = screen_center + rotated * (size * 0.5);

    var out: BillboardVertexOut;
    out.clip_position = screen_to_ndc(screen_pos);
    out.local_pos = vin.local_pos;
    out.color = iin.color;
    out.kind = iin.kind;
    return out;
}

/// Unsigned distance from `pt` to the segment `a`-`b` - used to build an
/// approximate signed distance to the triangle's boundary for
/// antialiasing (exact enough at the ~1px scale `aa_alpha` softens over).
fn distance_to_segment(pt: vec2<f32>, a: vec2<f32>, b: vec2<f32>) -> f32 {
    let ab = b - a;
    let t = clamp(dot(pt - a, ab) / max(dot(ab, ab), 1e-6), 0.0, 1.0);
    return length(pt - (a + ab * t));
}

// Standard sign-of-cross-product point-in-triangle test.
fn triangle_sign(p1: vec2<f32>, p2: vec2<f32>, p3: vec2<f32>) -> f32 {
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}

fn point_in_triangle(pt: vec2<f32>, v1: vec2<f32>, v2: vec2<f32>, v3: vec2<f32>) -> bool {
    let d1 = triangle_sign(pt, v1, v2);
    let d2 = triangle_sign(pt, v2, v3);
    let d3 = triangle_sign(pt, v3, v1);
    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    return !(has_neg && has_pos);
}

/// Signed approximate distance to the triangle's boundary: negative
/// inside, positive outside, magnitude the (approximate, for interior
/// points) distance to the nearest edge - enough for `aa_alpha` to soften
/// against, not a mathematically exact SDF.
fn triangle_sdf(pt: vec2<f32>, v1: vec2<f32>, v2: vec2<f32>, v3: vec2<f32>) -> f32 {
    let d = min(min(distance_to_segment(pt, v1, v2), distance_to_segment(pt, v2, v3)), distance_to_segment(pt, v3, v1));
    return select(d, -d, point_in_triangle(pt, v1, v2, v3));
}

@fragment
fn fs_billboard(in: BillboardVertexOut) -> @location(0) vec4<f32> {
    var alpha = 1.0;
    if (in.kind == 0u) {
        // circle
        alpha = aa_alpha(length(in.local_pos) - 1.0);
    } else if (in.kind == 2u || in.kind == 3u) {
        // Fixed "points along +Y" triangle in local space (reservoir
        // symbol or direction arrow) - the vertex shader's rotation is
        // what actually points it along `direction` on screen; this test
        // always runs in the pre-rotation frame.
        let tip = vec2<f32>(0.0, 1.0);
        let base_l = vec2<f32>(-0.9, -0.6);
        let base_r = vec2<f32>(0.9, -0.6);
        alpha = aa_alpha(triangle_sdf(in.local_pos, tip, base_l, base_r));
    }
    // kind == 1 (square): the whole quad is the shape, no discard/AA needed.
    if (alpha <= 0.0) {
        discard;
    }
    return vec4<f32>(in.color.rgb, in.color.a * alpha);
}

// --- links: thick lines ---
//
// `local_pos.x` in `[0, 1]` runs along the line (0 = start, 1 = end);
// `local_pos.y` in `[-1, 1]` runs across it (perpendicular offset in
// half-widths) - the standard "instanced quad as a thick line" technique.

struct LinkVertexIn {
    @location(0) local_pos: vec2<f32>,
};

struct LinkInstanceIn {
    @location(1) start: vec2<f32>,  // world space (origin-relative - see `Camera`)
    @location(2) end: vec2<f32>,    // world space (origin-relative - see `Camera`)
    @location(3) width: f32,        // pixels
    @location(4) color: vec4<f32>,
};

struct LinkVertexOut {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) across: f32, // local_pos.y carried through for the fragment shader's edge AA
};

@vertex
fn vs_link(vin: LinkVertexIn, iin: LinkInstanceIn) -> LinkVertexOut {
    let start_screen = world_to_screen(iin.start);
    let end_screen = world_to_screen(iin.end);
    let delta = end_screen - start_screen;
    let len = length(delta);
    let unit_dir = select(vec2<f32>(1.0, 0.0), delta / max(len, 0.0001), len > 0.0001);
    let normal = vec2<f32>(-unit_dir.y, unit_dir.x);

    let along = start_screen + unit_dir * (vin.local_pos.x * len);
    let screen_pos = along + normal * (vin.local_pos.y * iin.width * 0.5);

    var out: LinkVertexOut;
    out.clip_position = screen_to_ndc(screen_pos);
    out.color = iin.color;
    out.across = vin.local_pos.y;
    return out;
}

@fragment
fn fs_link(in: LinkVertexOut) -> @location(0) vec4<f32> {
    let alpha = aa_alpha(abs(in.across) - 1.0);
    if (alpha <= 0.0) {
        discard;
    }
    return vec4<f32>(in.color.rgb, in.color.a * alpha);
}
