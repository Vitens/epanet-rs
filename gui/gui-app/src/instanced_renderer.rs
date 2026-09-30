//! True GPU-instanced rendering for the network's nodes and links: one
//! small per-shape vertex buffer (a unit quad), one instance buffer per
//! shape *type* holding every node/link's data, and exactly one draw call
//! per shape type - "render the whole network in one pass," the way a
//! game engine renders a large sprite/particle system, rather than
//! Galileo's `FeatureLayer` (CPU-tessellate each feature once, batch into
//! ~256KB GPU buffers, one draw call *per batch*). See `network_shader.wgsl`
//! for the actual vertex/fragment shaders.
//!
//! Node symbols and direction-arrow triangles share one "billboard"
//! instance format and pipeline (both are fixed-pixel-size shapes,
//! oriented per-instance) and live in the *same* instance buffer - arrows
//! packed right after nodes (`billboard[0..node_count]` = nodes,
//! `billboard[node_count..node_count + link_count]` = one arrow per link,
//! `arrow_index == link_index`). Links (thick lines) get their own
//! instance format/pipeline, since a line's geometry (screen-space width
//! along an arbitrary-length segment) doesn't fit the same quad shape as
//! a billboard.
//!
//! Sync strategy is the same as the Galileo version it replaced (see the
//! git history of `network_layers.rs`), just targeting flat GPU buffers
//! indexed by the network's own `usize` node/link index directly instead
//! of Galileo's opaque `FeatureId` - a real simplification this rewrite
//! gets for free, since a GPU instance buffer *is* just an array:
//! - `topology_version` changed -> full rebuild **on a background
//!   thread** (constructing the `Vec<BillboardInstance>`/`Vec<LinkInstance>`
//!   is the expensive part, same as the Galileo version's
//!   `FeatureLayer::add` loop; only plain data leaves the thread, the
//!   actual `wgpu::Buffer` is created back on the UI thread).
//! - `properties_version` changed -> `queue.write_buffer` at exactly the
//!   changed indices' byte offsets (`network.updated_nodes`/`updated_links`).
//! - Selection changed -> same, diffed against the previous frame.
//! - Results changed -> patch every node/link's color field once.

use std::sync::mpsc::{Receiver, TryRecvError, channel};
use std::thread;

use bytemuck::{Pod, Zeroable};
use eframe::egui;
use wgpu::util::DeviceExt;

use epanet_rs::model::link::LinkType;
use epanet_rs::model::network::Network;
use epanet_rs::model::node::NodeType;

use gui_core::{AppState, Selection};

const SHADER_SOURCE: &str = include_str!("network_shader.wgsl");

const NODE_RADIUS_PX: f32 = 4.5;
const NODE_RADIUS_SELECTED_PX: f32 = 6.0;
const LINK_WIDTH_PX: f32 = 1.6;
const LINK_WIDTH_SELECTED_PX: f32 = 3.0;
const ARROW_SIZE_PX: f32 = 11.0;

mod style {
    pub const PIPE: [f32; 4] = [140.0 / 255.0, 170.0 / 255.0, 200.0 / 255.0, 1.0];
    pub const PUMP: [f32; 4] = [230.0 / 255.0, 90.0 / 255.0, 90.0 / 255.0, 1.0];
    pub const VALVE: [f32; 4] = [230.0 / 255.0, 200.0 / 255.0, 80.0 / 255.0, 1.0];

    pub const JUNCTION: [f32; 4] = [120.0 / 255.0, 190.0 / 255.0, 250.0 / 255.0, 1.0];
    pub const TANK: [f32; 4] = [250.0 / 255.0, 170.0 / 255.0, 80.0 / 255.0, 1.0];
    pub const RESERVOIR: [f32; 4] = [120.0 / 255.0, 220.0 / 255.0, 140.0 / 255.0, 1.0];

    // Not white: the canvas background is light (an empty Galileo
    // backdrop was briefly, accidentally overwriting it with white - see
    // `canvas.rs`'s basemap-backdrop comment - but even fixed, a light
    // background is the more typical case), so a pure-white highlight
    // risks being genuinely invisible against it. Orange reads clearly as
    // "selected" against both light and dark backgrounds.
    pub const SELECTED: [f32; 4] = [1.0, 0.55, 0.1, 1.0];
    pub const ARROW: [f32; 4] = [1.0, 1.0, 1.0, 0.8];
}

/// Shape discriminant for the "billboard" pipeline - matches
/// `network_shader.wgsl`'s `kind` values exactly.
#[repr(u32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BillboardKind {
    Circle = 0,
    Square = 1,
    /// Reservoir symbol - a fixed-pixel-size triangle, always drawn
    /// regardless of anything about the network's current zoom (unlike
    /// `Arrow`, below).
    Triangle = 2,
    /// A link's direction-arrow triangle - visually the same shape as
    /// `Triangle`, but a distinct `kind` so the vertex shader can hide it
    /// when its link is short on screen (`world_length`, below) without
    /// also hiding actual reservoir symbols, which have no comparable
    /// "shrink when zoomed out" concept.
    Arrow = 3,
}

/// One node symbol or one direction-arrow triangle - see this module's
/// doc comment for why both share this format/buffer/pipeline.
///
/// `center`/`direction` are relative to `InstancedRenderer`'s floating
/// origin, not raw model/INP coordinates - see that struct's doc comment
/// for why (`f32` precision at real-world coordinate scales).
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct BillboardInstance {
    center: [f32; 2],
    direction: [f32; 2],
    size: f32,
    color: [f32; 4],
    kind: u32,
    /// `Arrow` only: the link's full route length in world units, so the
    /// vertex shader can compute its current on-screen length
    /// (`world_length * camera.zoom`) and hide the arrow if its link is
    /// too short on screen to be worth drawing - see
    /// `network_shader.wgsl`'s `MIN_ARROW_SCREEN_LENGTH_PX`. Unused
    /// (`0.0`) for every other kind.
    world_length: f32,
}

impl BillboardInstance {
    const UP: [f32; 2] = [0.0, 1.0];

    fn node(origin: (f64, f64), center: (f64, f64), kind: BillboardKind, size: f32, color: [f32; 4]) -> Self {
        Self {
            center: [(center.0 - origin.0) as f32, (center.1 - origin.1) as f32],
            direction: Self::UP,
            size,
            color,
            kind: kind as u32,
            world_length: 0.0,
        }
    }

    fn arrow(
        origin: (f64, f64),
        center: (f64, f64),
        direction: (f32, f32),
        world_length: f32,
        color: [f32; 4],
    ) -> Self {
        Self {
            center: [(center.0 - origin.0) as f32, (center.1 - origin.1) as f32],
            direction: [direction.0, direction.1],
            size: ARROW_SIZE_PX,
            color,
            kind: BillboardKind::Arrow as u32,
            world_length,
        }
    }

    /// A degenerate (zero-size, therefore invisible) instance - used to
    /// fill an arrow slot for a link whose start/end can't be resolved
    /// (missing coordinates), keeping `arrow_index == link_index`
    /// alignment without a secondary "has no arrow" side table.
    fn hidden() -> Self {
        Self {
            center: [f32::NAN, f32::NAN],
            direction: Self::UP,
            size: 0.0,
            color: [0.0, 0.0, 0.0, 0.0],
            kind: BillboardKind::Circle as u32,
            world_length: 0.0,
        }
    }
}

/// `start`/`end` are relative to `InstancedRenderer`'s floating origin,
/// not raw model/INP coordinates - see that struct's doc comment.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct LinkInstance {
    start: [f32; 2],
    end: [f32; 2],
    width: f32,
    color: [f32; 4],
    _pad: [f32; 1],
}

impl LinkInstance {
    fn new(origin: (f64, f64), start: (f64, f64), end: (f64, f64), width: f32, color: [f32; 4]) -> Self {
        Self {
            start: [(start.0 - origin.0) as f32, (start.1 - origin.1) as f32],
            end: [(end.0 - origin.0) as f32, (end.1 - origin.1) as f32],
            width,
            color,
            _pad: [0.0],
        }
    }

    fn hidden() -> Self {
        Self::new((0.0, 0.0), (f64::NAN, f64::NAN), (f64::NAN, f64::NAN), 0.0, [0.0; 4])
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct CameraUniform {
    center: [f32; 2],
    zoom: f32,
    _pad0: f32,
    viewport: [f32; 2],
    _pad1: [f32; 2],
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct QuadVertex {
    local_pos: [f32; 2],
}

const QUAD_INDICES: [u16; 6] = [0, 1, 2, 0, 2, 3];
/// Billboard quad: a centered unit square, `[-1, 1]` on both axes.
const BILLBOARD_QUAD: [QuadVertex; 4] = [
    QuadVertex { local_pos: [-1.0, -1.0] },
    QuadVertex { local_pos: [1.0, -1.0] },
    QuadVertex { local_pos: [1.0, 1.0] },
    QuadVertex { local_pos: [-1.0, 1.0] },
];
/// Link quad: `x` runs `[0, 1]` along the line, `y` runs `[-1, 1]` across
/// it (half-widths) - see `network_shader.wgsl::vs_link`.
const LINK_QUAD: [QuadVertex; 4] = [
    QuadVertex { local_pos: [0.0, -1.0] },
    QuadVertex { local_pos: [1.0, -1.0] },
    QuadVertex { local_pos: [1.0, 1.0] },
    QuadVertex { local_pos: [0.0, 1.0] },
];

fn instance_vertex_buffer_layout<'a>(
    array_stride: u64,
    attributes: &'a [wgpu::VertexAttribute],
) -> wgpu::VertexBufferLayout<'a> {
    wgpu::VertexBufferLayout {
        array_stride,
        step_mode: wgpu::VertexStepMode::Instance,
        attributes,
    }
}

/// An instance buffer that gets recreated wholesale on a topology rebuild
/// and patched in place (via `queue.write_buffer` at one instance's byte
/// offset) for property/selection/result updates - no batching, no
/// "updating one invalidates its neighbors" concern at all, unlike
/// Galileo's bundle system: a GPU buffer is just an array, and the
/// network's own index *is* the array index.
struct InstanceBuffer {
    buffer: Option<wgpu::Buffer>,
    count: usize,
}

impl InstanceBuffer {
    fn empty() -> Self {
        Self { buffer: None, count: 0 }
    }

    fn upload<T: Pod>(&mut self, device: &wgpu::Device, label: &str, data: &[T]) {
        self.count = data.len();
        self.buffer = if data.is_empty() {
            None
        } else {
            Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: bytemuck::cast_slice(data),
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            }))
        };
    }

    fn patch<T: Pod>(&self, queue: &wgpu::Queue, index: usize, value: &T) {
        if let Some(buffer) = &self.buffer {
            let offset = (index * std::mem::size_of::<T>()) as u64;
            queue.write_buffer(buffer, offset, bytemuck::bytes_of(value));
        }
    }
}

struct PendingRebuild {
    target_topology_version: u32,
    target_properties_version: u32,
    target_selected_nodes: Vec<String>,
    target_selected_links: Vec<String>,
    receiver: Receiver<RebuildOutput>,
}

struct RebuildOutput {
    nodes: Vec<BillboardInstance>,
    links: Vec<LinkInstance>,
    link_segment_offsets: Vec<usize>,
    arrows: Vec<BillboardInstance>,
    origin: (f64, f64),
}

pub struct InstancedRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,

    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,

    billboard_pipeline: wgpu::RenderPipeline,
    billboard_vertex_buffer: wgpu::Buffer,
    link_pipeline: wgpu::RenderPipeline,
    link_vertex_buffer: wgpu::Buffer,
    quad_index_buffer: wgpu::Buffer,

    /// Nodes (`[0..node_count]`) and one arrow per link
    /// (`[node_count..node_count + link_count]`) - see the module doc
    /// comment for why they share a buffer/pipeline.
    billboards: InstanceBuffer,
    node_count: usize,
    links: InstanceBuffer,
    /// `network.links[i]`'s segments live at
    /// `links[link_segment_offsets[i]..link_segment_offsets[i + 1]]` - see
    /// `build_all`'s doc comment for why a link can span more than one
    /// `LinkInstance` (multi-vertex pipes).
    link_segment_offsets: Vec<usize>,

    /// A fixed reference point (some node/link endpoint's coordinates,
    /// picked once per rebuild - see `build_all`), subtracted from every
    /// world-space value *before* it's cast to the `f32` the GPU buffers
    /// and shaders use. Real networks use real-world coordinate systems
    /// (Dutch RD grid puts values in the hundreds of thousands of meters);
    /// `f32` has ~7 significant decimal digits, so storing those coordinates
    /// directly leaves only ~1 digit of sub-meter precision - visibly
    /// wrong/jittery node positions and pipe angles on a dense real network.
    /// Subtracting a nearby origin first keeps every stored/uniform value
    /// small (typically a few thousand at most), where `f32` has ample
    /// precision. `set_camera` subtracts the same origin from the camera
    /// center for the same reason, and the shader never sees raw model
    /// coordinates at all.
    origin: (f64, f64),

    topology_version: Option<u32>,
    properties_version: u32,
    selected_nodes: Vec<String>,
    selected_links: Vec<String>,
    results_signature: Option<(usize, usize)>,
    pending: Option<PendingRebuild>,
}

impl InstancedRenderer {
    pub fn new(render_state: &egui_wgpu::RenderState, state: &AppState) -> Self {
        let device = render_state.device.clone();
        let queue = render_state.queue.clone();

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("network_shader"),
            source: wgpu::ShaderSource::Wgsl(SHADER_SOURCE.into()),
        });

        let camera_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("network_camera"),
            size: std::mem::size_of::<CameraUniform>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let camera_bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("network_camera_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });
        let camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("network_camera_bind_group"),
            layout: &camera_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::Buffer(camera_buffer.as_entire_buffer_binding()),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("network_pipeline_layout"),
            bind_group_layouts: &[Some(&camera_bind_group_layout)],
            immediate_size: 0,
        });

        let target_format = render_state.target_format;
        let color_target = Some(wgpu::ColorTargetState {
            format: target_format,
            blend: Some(wgpu::BlendState::ALPHA_BLENDING),
            write_mask: wgpu::ColorWrites::ALL,
        });

        let vertex_attr = wgpu::vertex_attr_array![0 => Float32x2];
        let billboard_instance_attrs = wgpu::vertex_attr_array![
            1 => Float32x2, // center
            2 => Float32x2, // direction
            3 => Float32,   // size
            4 => Float32x4, // color
            5 => Uint32,    // kind
            6 => Float32,   // world_length
        ];
        let link_instance_attrs = wgpu::vertex_attr_array![
            1 => Float32x2, // start
            2 => Float32x2, // end
            3 => Float32,   // width
            4 => Float32x4, // color
        ];

        let billboard_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("network_billboard_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_billboard"),
                buffers: &[
                    wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<QuadVertex>() as u64,
                        step_mode: wgpu::VertexStepMode::Vertex,
                        attributes: &vertex_attr,
                    },
                    instance_vertex_buffer_layout(
                        std::mem::size_of::<BillboardInstance>() as u64,
                        &billboard_instance_attrs,
                    ),
                ],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_billboard"),
                targets: std::slice::from_ref(&color_target),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview_mask: None,
            cache: None,
        });

        let link_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("network_link_pipeline"),
            layout: Some(&pipeline_layout),
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_link"),
                buffers: &[
                    wgpu::VertexBufferLayout {
                        array_stride: std::mem::size_of::<QuadVertex>() as u64,
                        step_mode: wgpu::VertexStepMode::Vertex,
                        attributes: &vertex_attr,
                    },
                    instance_vertex_buffer_layout(
                        std::mem::size_of::<LinkInstance>() as u64,
                        &link_instance_attrs,
                    ),
                ],
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
                entry_point: Some("fs_link"),
                targets: std::slice::from_ref(&color_target),
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            }),
            primitive: wgpu::PrimitiveState::default(),
            depth_stencil: None,
            multisample: wgpu::MultisampleState::default(),
            multiview_mask: None,
            cache: None,
        });

        let billboard_vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("network_billboard_quad"),
            contents: bytemuck::cast_slice(&BILLBOARD_QUAD),
            usage: wgpu::BufferUsages::VERTEX,
        });
        let link_vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("network_link_quad"),
            contents: bytemuck::cast_slice(&LINK_QUAD),
            usage: wgpu::BufferUsages::VERTEX,
        });
        let quad_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("network_quad_indices"),
            contents: bytemuck::cast_slice(&QUAD_INDICES),
            usage: wgpu::BufferUsages::INDEX,
        });

        let mut renderer = Self {
            device,
            queue,
            camera_buffer,
            camera_bind_group,
            billboard_pipeline,
            billboard_vertex_buffer,
            link_pipeline,
            link_vertex_buffer,
            quad_index_buffer,
            billboards: InstanceBuffer::empty(),
            node_count: 0,
            links: InstanceBuffer::empty(),
            link_segment_offsets: Vec::new(),
            origin: (0.0, 0.0),
            topology_version: None,
            properties_version: 0,
            selected_nodes: state.selection.nodes.clone(),
            selected_links: state.selection.links.clone(),
            results_signature: None,
            pending: None,
        };
        renderer.spawn_rebuild(state);
        renderer
    }

    pub fn set_camera(&self, center: (f64, f64), zoom: f64, viewport: (f32, f32)) {
        let relative = (center.0 - self.origin.0, center.1 - self.origin.1);
        let uniform = CameraUniform {
            center: [relative.0 as f32, relative.1 as f32],
            zoom: zoom as f32,
            _pad0: 0.0,
            viewport: [viewport.0, viewport.1],
            _pad1: [0.0, 0.0],
        };
        self.queue.write_buffer(&self.camera_buffer, 0, bytemuck::bytes_of(&uniform));
    }

    /// Mirrors any changes in `state` into the GPU instance buffers.
    /// Cheap (no-op) when nothing relevant changed since the last call.
    pub fn sync(&mut self, state: &AppState) {
        let network = state.network();

        if let Some(pending) = self.pending.take() {
            match pending.receiver.try_recv() {
                Ok(output) => {
                    let stale = pending.target_topology_version != network.topology_version;
                    if stale {
                        self.spawn_rebuild(state);
                    } else {
                        self.node_count = output.nodes.len();
                        let mut billboards = output.nodes;
                        billboards.extend(output.arrows);
                        self.billboards.upload(&self.device, "network_billboards", &billboards);
                        self.links.upload(&self.device, "network_links", &output.links);
                        self.link_segment_offsets = output.link_segment_offsets;
                        self.origin = output.origin;

                        self.topology_version = Some(network.topology_version);
                        self.properties_version = pending.target_properties_version;
                        self.selected_nodes = pending.target_selected_nodes;
                        self.selected_links = pending.target_selected_links;
                        self.results_signature = None;
                    }
                }
                Err(TryRecvError::Empty) => {
                    self.pending = Some(pending);
                    return;
                }
                Err(TryRecvError::Disconnected) => {}
            }
        }

        if self.topology_version != Some(network.topology_version) {
            self.spawn_rebuild(state);
            return;
        }

        if network.properties_version != self.properties_version {
            self.sync_properties(state);
            self.properties_version = network.properties_version;
        }

        self.sync_selection(state);

        let results_signature = results_signature(state);
        if results_signature != self.results_signature {
            self.sync_results(state);
            self.results_signature = results_signature;
        }
    }

    fn spawn_rebuild(&mut self, state: &AppState) {
        let network = state.network().clone();
        let selection = state.selection.clone();
        let target_topology_version = network.topology_version;
        let target_properties_version = network.properties_version;
        let target_selected_nodes = selection.nodes.clone();
        let target_selected_links = selection.links.clone();

        let (tx, rx) = channel();
        thread::spawn(move || {
            let (nodes, links, link_segment_offsets, arrows, origin) = build_all(&network, &selection, None);
            let _ = tx.send(RebuildOutput { nodes, links, link_segment_offsets, arrows, origin });
        });

        self.pending = Some(PendingRebuild {
            target_topology_version,
            target_properties_version,
            target_selected_nodes,
            target_selected_links,
            receiver: rx,
        });
    }

    /// Patches every segment of `link_index`'s route (see `link_points`/
    /// `build_all`) in place. If the link's vertex count changed since the
    /// last rebuild (editing a pipe's display route, rather than moving
    /// an endpoint or changing a color-relevant property), the number of
    /// segments needed may no longer match the slots
    /// `link_segment_offsets` allocated for it at that rebuild - extra
    /// segments are dropped and missing ones leave stale geometry in
    /// their slot until the next full rebuild, rather than growing the
    /// buffer on the spot. A real, deliberate limitation: vertex-count
    /// edits aren't currently exposed as their own action anyway (there's
    /// no `Action::SetLinkVertices`), so this path isn't reachable from
    /// the UI today - flagged here so it isn't a silent surprise if that
    /// changes.
    fn patch_link(&self, network: &Network, link_index: usize, selected: bool, result_color: Option<[f32; 4]>) {
        let Some(link) = network.links.get(link_index) else {
            return;
        };
        let (Some(&start_offset), Some(&end_offset)) = (
            self.link_segment_offsets.get(link_index),
            self.link_segment_offsets.get(link_index + 1),
        ) else {
            return;
        };

        match link_points(network, link) {
            Some(points) => {
                for (i, instance) in link_segment_instances(self.origin, link, &points, selected, result_color)
                    .into_iter()
                    .enumerate()
                {
                    let slot = start_offset + i;
                    if slot >= end_offset {
                        break;
                    }
                    self.links.patch(&self.queue, slot, &instance);
                }
            }
            None => {
                for slot in start_offset..end_offset {
                    self.links.patch(&self.queue, slot, &LinkInstance::hidden());
                }
            }
        }
    }

    fn patch_arrow(&self, network: &Network, link_index: usize, reversed: bool) {
        let Some(link) = network.links.get(link_index) else {
            return;
        };
        let arrow = match link_points(network, link) {
            Some(points) => arrow_instance(self.origin, &points, reversed),
            None => BillboardInstance::hidden(),
        };
        self.billboards.patch(&self.queue, self.node_count + link_index, &arrow);
    }

    fn sync_properties(&mut self, state: &AppState) {
        let network = state.network();

        for &node_index in &network.updated_nodes {
            let Some(node) = network.nodes.get(node_index) else {
                continue;
            };
            let instance = node_instance(self.origin, node, state.selection.contains_node(&node.id));
            self.billboards.patch(&self.queue, node_index, &instance);
        }

        for &link_index in &network.updated_links {
            let Some(link) = network.links.get(link_index) else {
                continue;
            };
            let selected = state.selection.contains_link(&link.id);
            self.patch_link(network, link_index, selected, None);
            self.patch_arrow(network, link_index, false);
        }
    }

    fn sync_selection(&mut self, state: &AppState) {
        let network = state.network();

        if state.selection.nodes != self.selected_nodes {
            for id in self.selected_nodes.iter().chain(state.selection.nodes.iter()) {
                let Some(&node_index) = network.node_map.get(id.as_str()) else {
                    continue;
                };
                let Some(node) = network.nodes.get(node_index) else {
                    continue;
                };
                let instance = node_instance(self.origin, node, state.selection.contains_node(id));
                self.billboards.patch(&self.queue, node_index, &instance);
            }
            self.selected_nodes = state.selection.nodes.clone();
        }

        if state.selection.links != self.selected_links {
            for id in self.selected_links.iter().chain(state.selection.links.iter()) {
                let Some(&link_index) = network.link_map.get(id.as_str()) else {
                    continue;
                };
                let selected = state.selection.contains_link(id);
                self.patch_link(network, link_index, selected, None);
            }
            self.selected_links = state.selection.links.clone();
        }
    }

    fn sync_results(&mut self, state: &AppState) {
        let network = state.network();
        let overlay = super::canvas::ResultOverlay::new(state);

        for (node_index, node) in network.nodes.iter().enumerate() {
            let color = overlay.node_color_rgb(node_index).map(rgb_to_f32);
            let mut instance = node_instance(self.origin, node, state.selection.contains_node(&node.id));
            if let Some(color) = color {
                instance.color = color;
            }
            self.billboards.patch(&self.queue, node_index, &instance);
        }

        for (link_index, link) in network.links.iter().enumerate() {
            let selected = state.selection.contains_link(&link.id);
            let color = overlay.link_color_rgb(link_index).map(rgb_to_f32);
            self.patch_link(network, link_index, selected, color);

            let reversed = overlay.flow_reversed(link_index);
            self.patch_arrow(network, link_index, reversed);
        }
    }

    /// Builds an `egui::Shape` that issues this frame's two draw calls
    /// (billboards, links) - one `egui_wgpu::Callback` each, painted in
    /// `rect`. Cheap to construct: it just clones a handful of wgpu
    /// handles (all internally `Arc`-based) into small per-frame structs.
    pub fn paint_shape(&self, rect: egui::Rect) -> Vec<egui::Shape> {
        let mut shapes = Vec::new();

        if let (Some(buffer), true) = (&self.links.buffer, self.links.count > 0) {
            shapes.push(
                egui_wgpu::Callback::new_paint_callback(
                    rect,
                    LinkDraw {
                        pipeline: self.link_pipeline.clone(),
                        camera_bind_group: self.camera_bind_group.clone(),
                        vertex_buffer: self.link_vertex_buffer.clone(),
                        index_buffer: self.quad_index_buffer.clone(),
                        instance_buffer: buffer.clone(),
                        instance_count: self.links.count as u32,
                    },
                )
                .into(),
            );
        }

        if let (Some(buffer), true) = (&self.billboards.buffer, self.billboards.count > 0) {
            shapes.push(
                egui_wgpu::Callback::new_paint_callback(
                    rect,
                    BillboardDraw {
                        pipeline: self.billboard_pipeline.clone(),
                        camera_bind_group: self.camera_bind_group.clone(),
                        vertex_buffer: self.billboard_vertex_buffer.clone(),
                        index_buffer: self.quad_index_buffer.clone(),
                        instance_buffer: buffer.clone(),
                        instance_count: self.billboards.count as u32,
                    },
                )
                .into(),
            );
        }

        shapes
    }
}

struct BillboardDraw {
    pipeline: wgpu::RenderPipeline,
    camera_bind_group: wgpu::BindGroup,
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    instance_buffer: wgpu::Buffer,
    instance_count: u32,
}

impl egui_wgpu::CallbackTrait for BillboardDraw {
    fn paint(
        &self,
        _info: egui::epaint::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        _resources: &egui_wgpu::CallbackResources,
    ) {
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
        render_pass.set_vertex_buffer(1, self.instance_buffer.slice(..));
        render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
        render_pass.draw_indexed(0..QUAD_INDICES.len() as u32, 0, 0..self.instance_count);
    }
}

struct LinkDraw {
    pipeline: wgpu::RenderPipeline,
    camera_bind_group: wgpu::BindGroup,
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    instance_buffer: wgpu::Buffer,
    instance_count: u32,
}

impl egui_wgpu::CallbackTrait for LinkDraw {
    fn paint(
        &self,
        _info: egui::epaint::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        _resources: &egui_wgpu::CallbackResources,
    ) {
        render_pass.set_pipeline(&self.pipeline);
        render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.vertex_buffer.slice(..));
        render_pass.set_vertex_buffer(1, self.instance_buffer.slice(..));
        render_pass.set_index_buffer(self.index_buffer.slice(..), wgpu::IndexFormat::Uint16);
        render_pass.draw_indexed(0..QUAD_INDICES.len() as u32, 0, 0..self.instance_count);
    }
}

fn rgb_to_f32((r, g, b): (u8, u8, u8)) -> [f32; 4] {
    [r as f32 / 255.0, g as f32 / 255.0, b as f32 / 255.0, 1.0]
}

fn node_instance(origin: (f64, f64), node: &epanet_rs::model::node::Node, selected: bool) -> BillboardInstance {
    let Some(coordinates) = node.coordinates else {
        return BillboardInstance::hidden();
    };
    let kind = match node.node_type {
        NodeType::Junction(_) => BillboardKind::Circle,
        NodeType::Tank(_) => BillboardKind::Square,
        NodeType::Reservoir(_) => BillboardKind::Triangle,
    };
    let base_color = match node.node_type {
        NodeType::Junction(_) => style::JUNCTION,
        NodeType::Tank(_) => style::TANK,
        NodeType::Reservoir(_) => style::RESERVOIR,
    };
    let (size, color) = if selected {
        (NODE_RADIUS_SELECTED_PX * 2.0, style::SELECTED)
    } else {
        (NODE_RADIUS_PX * 2.0, base_color)
    };
    BillboardInstance::node(origin, coordinates, kind, size, color)
}

fn link_endpoints(network: &Network, link: &epanet_rs::model::link::Link) -> Option<((f64, f64), (f64, f64))> {
    let start = network.nodes.get(link.start_node)?.coordinates?;
    let end = network.nodes.get(link.end_node)?.coordinates?;
    Some((start, end))
}

/// The link's full displayed route - `[start, ...link.vertices, end]` -
/// matching exactly what the earlier Galileo version's `build_link_line`
/// drew (see the module doc: this replaces the instanced renderer's
/// earlier single-straight-chord simplification). `None` if either
/// endpoint node has no coordinates.
fn link_points(network: &Network, link: &epanet_rs::model::link::Link) -> Option<Vec<(f64, f64)>> {
    let (start, end) = link_endpoints(network, link)?;
    let mut points = Vec::with_capacity(2 + link.vertices.as_ref().map_or(0, Vec::len));
    points.push(start);
    if let Some(vertices) = &link.vertices {
        points.extend(vertices.iter().copied());
    }
    points.push(end);
    Some(points)
}

fn link_color(link: &epanet_rs::model::link::Link, selected: bool, result_color: Option<[f32; 4]>) -> (f32, [f32; 4]) {
    let base_color = result_color.unwrap_or(match link.link_type {
        LinkType::Pipe(_) => style::PIPE,
        LinkType::Pump(_) => style::PUMP,
        LinkType::Valve(_) => style::VALVE,
    });
    if selected {
        (LINK_WIDTH_SELECTED_PX, style::SELECTED)
    } else {
        (LINK_WIDTH_PX, base_color)
    }
}

/// One `LinkInstance` per segment of `points` (a route of 2+ points) - see
/// `link_points`. A straight two-point link produces exactly one
/// instance, same as before; a multi-vertex one produces one per segment,
/// all sharing the same width/color.
fn link_segment_instances(
    origin: (f64, f64),
    link: &epanet_rs::model::link::Link,
    points: &[(f64, f64)],
    selected: bool,
    result_color: Option<[f32; 4]>,
) -> Vec<LinkInstance> {
    let (width, color) = link_color(link, selected, result_color);
    points
        .windows(2)
        .map(|w| LinkInstance::new(origin, w[0], w[1], width, color))
        .collect()
}

/// Walks `points` (a route of 2+ points, see `link_points`) to its
/// halfway point by total length, tracking the local tangent of whichever
/// segment it falls on - so the direction arrow lands correctly on a
/// multi-vertex pipe's actual midpoint, not the straight chord between its
/// two endpoints. Ported from the earlier Galileo version's
/// `direction_arrow` (see git history), which had the same requirement.
/// Also returns the route's total length (`world_length`, for the
/// shader's screen-space visibility cull - see
/// `BillboardInstance::world_length`).
fn arrow_instance(origin: (f64, f64), points: &[(f64, f64)], reversed: bool) -> BillboardInstance {
    if points.len() < 2 {
        return BillboardInstance::hidden();
    }
    let seg_len = |a: (f64, f64), b: (f64, f64)| {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        (dx * dx + dy * dy).sqrt()
    };
    let total: f64 = points.windows(2).map(|w| seg_len(w[0], w[1])).sum();
    if total < 1e-9 {
        return BillboardInstance::hidden();
    }

    let mut remaining = total * 0.5;
    let mut mid = points[0];
    let mut tangent = (1.0, 0.0);
    for w in points.windows(2) {
        let length = seg_len(w[0], w[1]);
        if length < 1e-12 {
            continue;
        }
        let (dx, dy) = ((w[1].0 - w[0].0) / length, (w[1].1 - w[0].1) / length);
        if remaining <= length {
            mid = (w[0].0 + dx * remaining, w[0].1 + dy * remaining);
            tangent = (dx, dy);
            break;
        }
        remaining -= length;
        mid = w[1];
        tangent = (dx, dy);
    }

    let (mut ux, mut uy) = (tangent.0 as f32, tangent.1 as f32);
    if reversed {
        ux = -ux;
        uy = -uy;
    }
    BillboardInstance::arrow(origin, mid, (ux, uy), total as f32, style::ARROW)
}

fn results_signature(state: &AppState) -> Option<(usize, usize)> {
    state
        .sim_results
        .as_ref()
        .map(|_| (state.report_step, state.network().topology_version as usize))
}

/// Picks a floating origin for a rebuild - see `InstancedRenderer::origin`'s
/// doc comment. Just the first node with coordinates (falling back to the
/// first link's start node if no node has one, then `(0.0, 0.0)` for a
/// fully coordinate-less network): it doesn't need to be the network's
/// true centroid or anything more principled than "some point actually
/// near the network's real coordinates," since its only job is keeping
/// every *other* coordinate's distance from it small enough for `f32` to
/// represent precisely.
fn pick_origin(network: &Network) -> (f64, f64) {
    network
        .nodes
        .iter()
        .find_map(|node| node.coordinates)
        .unwrap_or((0.0, 0.0))
}

/// Builds the full instance data for the current network - runs on
/// `spawn_rebuild`'s background thread. Also returns `link_segment_offsets`
/// (length `network.links.len() + 1`, cumulative): link `i`'s segments
/// live at `links[link_segment_offsets[i]..link_segment_offsets[i + 1]]`,
/// since a multi-vertex link occupies more than one `LinkInstance` (see
/// `link_points`) - every later per-link patch (`sync_properties`,
/// `sync_selection`, `sync_results`) needs this to find the right range.
/// Also returns the floating origin picked for this rebuild (see
/// `pick_origin`) - every instance here is already relative to it.
///
/// `_overlay` is always `None` here (results aren't threaded through the
/// background rebuild at all - `InstancedRenderer::sync` stamps the landed
/// buffers' `results_signature` back to `None`, so if results are already
/// active, the normal `sync_results` fast path picks them up immediately
/// afterward instead of this function needing to duplicate that logic).
fn build_all(
    network: &Network,
    selection: &Selection,
    _overlay: Option<()>,
) -> (Vec<BillboardInstance>, Vec<LinkInstance>, Vec<usize>, Vec<BillboardInstance>, (f64, f64)) {
    let origin = pick_origin(network);

    let nodes: Vec<BillboardInstance> = network
        .nodes
        .iter()
        .map(|node| node_instance(origin, node, selection.contains_node(&node.id)))
        .collect();

    let mut links = Vec::new();
    let mut arrows = Vec::with_capacity(network.links.len());
    let mut link_segment_offsets = Vec::with_capacity(network.links.len() + 1);
    for link in &network.links {
        link_segment_offsets.push(links.len());
        match link_points(network, link) {
            Some(points) => {
                let selected = selection.contains_link(&link.id);
                links.extend(link_segment_instances(origin, link, &points, selected, None));
                arrows.push(arrow_instance(origin, &points, false));
            }
            None => {
                links.push(LinkInstance::hidden());
                arrows.push(BillboardInstance::hidden());
            }
        }
    }
    link_segment_offsets.push(links.len());

    (nodes, links, link_segment_offsets, arrows, origin)
}
