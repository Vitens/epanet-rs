//! Uniform-grid spatial index over `Network` nodes and links.
//!
//! Rebuilt whenever `Network.topology_version` (or node/link coordinates)
//! change. Used for both frustum culling (in a renderer) and hit-testing
//! (click/hover/box-select), per the plan: the same index backs both.
//!
//! This is a dependency-free stand-in for an R-tree (`rstar`); a uniform grid
//! is simple, allocation-light to rebuild, and plenty fast for the node/link
//! counts this editor targets.

use epanet_rs::model::network::Network;

use crate::camera::Bounds;

type Cell = (i64, i64);

#[derive(Debug, Default, Clone)]
pub struct SpatialIndex {
    cell_size: f64,
    node_cells: hashbrown::HashMap<Cell, Vec<usize>>,
    link_cells: hashbrown::HashMap<Cell, Vec<usize>>,
    /// Cached node point per node index (only nodes with coordinates are indexed).
    node_points: Vec<Option<(f64, f64)>>,
    /// Cached polyline per link index (endpoints + vertices), in world space.
    link_lines: Vec<Option<Vec<(f64, f64)>>>,
}

fn cell_of(x: f64, y: f64, cell_size: f64) -> Cell {
    (
        (x / cell_size).floor() as i64,
        (y / cell_size).floor() as i64,
    )
}

impl SpatialIndex {
    /// Rebuild the index from scratch. `cell_size` should be roughly the
    /// scale of a typical link length; 100.0 is a reasonable default for
    /// most INP files (feet/meters-scale coordinates).
    pub fn build(network: &Network, cell_size: f64) -> Self {
        let cell_size = if cell_size > 0.0 { cell_size } else { 100.0 };
        let mut index = SpatialIndex {
            cell_size,
            node_cells: hashbrown::HashMap::new(),
            link_cells: hashbrown::HashMap::new(),
            node_points: Vec::with_capacity(network.nodes.len()),
            link_lines: Vec::with_capacity(network.links.len()),
        };

        for (i, node) in network.nodes.iter().enumerate() {
            index.node_points.push(node.coordinates);
            if let Some((x, y)) = node.coordinates {
                index
                    .node_cells
                    .entry(cell_of(x, y, cell_size))
                    .or_default()
                    .push(i);
            }
        }

        for (i, link) in network.links.iter().enumerate() {
            let start = network
                .nodes
                .get(link.start_node)
                .and_then(|n| n.coordinates);
            let end = network.nodes.get(link.end_node).and_then(|n| n.coordinates);
            let line = match (start, end) {
                (Some(s), Some(e)) => {
                    let mut pts = vec![s];
                    if let Some(vertices) = &link.vertices {
                        pts.extend(vertices.iter().copied());
                    }
                    pts.push(e);
                    Some(pts)
                }
                _ => None,
            };
            if let Some(pts) = &line {
                for cell in cells_touched_by_line(pts, cell_size) {
                    index.link_cells.entry(cell).or_default().push(i);
                }
            }
            index.link_lines.push(line);
        }

        index
    }

    pub fn node_point(&self, node_index: usize) -> Option<(f64, f64)> {
        self.node_points.get(node_index).copied().flatten()
    }

    pub fn link_line(&self, link_index: usize) -> Option<&[(f64, f64)]> {
        self.link_lines.get(link_index).and_then(|l| l.as_deref())
    }

    /// Find the closest node to `point` within `max_dist` (world units).
    pub fn nearest_node(&self, point: (f64, f64), max_dist: f64) -> Option<usize> {
        let (cx, cy) = cell_of(point.0, point.1, self.cell_size);
        let radius = (max_dist / self.cell_size).ceil() as i64 + 1;
        let mut best: Option<(usize, f64)> = None;
        for gx in (cx - radius)..=(cx + radius) {
            for gy in (cy - radius)..=(cy + radius) {
                if let Some(indices) = self.node_cells.get(&(gx, gy)) {
                    for &i in indices {
                        if let Some((x, y)) = self.node_points[i] {
                            let d = ((x - point.0).powi(2) + (y - point.1).powi(2)).sqrt();
                            if d <= max_dist && best.map(|(_, bd)| d < bd).unwrap_or(true) {
                                best = Some((i, d));
                            }
                        }
                    }
                }
            }
        }
        best.map(|(i, _)| i)
    }

    /// Find the closest link to `point` within `max_dist` (world units),
    /// measured as the distance from `point` to the link's polyline.
    pub fn nearest_link(&self, point: (f64, f64), max_dist: f64) -> Option<usize> {
        let (cx, cy) = cell_of(point.0, point.1, self.cell_size);
        let radius = (max_dist / self.cell_size).ceil() as i64 + 1;
        let mut seen: hashbrown::HashSet<usize> = hashbrown::HashSet::new();
        let mut best: Option<(usize, f64)> = None;
        for gx in (cx - radius)..=(cx + radius) {
            for gy in (cy - radius)..=(cy + radius) {
                if let Some(indices) = self.link_cells.get(&(gx, gy)) {
                    for &i in indices {
                        if !seen.insert(i) {
                            continue;
                        }
                        if let Some(line) = &self.link_lines[i] {
                            let d = distance_to_polyline(point, line);
                            if d <= max_dist && best.map(|(_, bd)| d < bd).unwrap_or(true) {
                                best = Some((i, d));
                            }
                        }
                    }
                }
            }
        }
        best.map(|(i, _)| i)
    }

    /// All node indices whose point lies within `bounds` (for box-select).
    pub fn nodes_in_bounds(&self, bounds: Bounds) -> Vec<usize> {
        self.node_points
            .iter()
            .enumerate()
            .filter_map(|(i, p)| {
                let (x, y) = (*p)?;
                bounds.contains_point(x, y).then_some(i)
            })
            .collect()
    }

    /// All node indices visible within `bounds`, for renderer culling.
    pub fn visible_nodes(&self, bounds: &Bounds) -> Vec<usize> {
        self.node_points
            .iter()
            .enumerate()
            .filter_map(|(i, p)| {
                let (x, y) = (*p)?;
                bounds.contains_point(x, y).then_some(i)
            })
            .collect()
    }

    /// All link indices whose bounding box intersects `bounds`, for renderer
    /// culling.
    pub fn visible_links(&self, bounds: &Bounds) -> Vec<usize> {
        self.link_lines
            .iter()
            .enumerate()
            .filter_map(|(i, line)| {
                let line = line.as_ref()?;
                let line_bounds = Bounds::from_points(line.iter().copied())?;
                bounds.intersects(&line_bounds).then_some(i)
            })
            .collect()
    }
}

fn cells_touched_by_line(points: &[(f64, f64)], cell_size: f64) -> Vec<Cell> {
    let mut cells = Vec::new();
    for window in points.windows(2) {
        let (x0, y0) = window[0];
        let (x1, y1) = window[1];
        // Coarse rasterization: step along the segment in cell-sized
        // increments. Fine for the moderate segment lengths typical of pipe
        // vertices; not a true Bresenham/DDA line but sufficient for
        // culling/hit-test bucketing purposes.
        let len = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
        let steps = ((len / cell_size).ceil() as usize).max(1);
        for step in 0..=steps {
            let t = step as f64 / steps as f64;
            let x = x0 + (x1 - x0) * t;
            let y = y0 + (y1 - y0) * t;
            cells.push(cell_of(x, y, cell_size));
        }
    }
    cells.sort_unstable();
    cells.dedup();
    cells
}

fn distance_to_polyline(point: (f64, f64), line: &[(f64, f64)]) -> f64 {
    line.windows(2)
        .map(|w| distance_to_segment(point, w[0], w[1]))
        .fold(f64::MAX, f64::min)
}

fn distance_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (px, py) = p;
    let (ax, ay) = a;
    let (bx, by) = b;
    let (dx, dy) = (bx - ax, by - ay);
    let len_sq = dx * dx + dy * dy;
    let t = if len_sq > 0.0 {
        (((px - ax) * dx + (py - ay) * dy) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let (cx, cy) = (ax + t * dx, ay + t * dy);
    ((px - cx).powi(2) + (py - cy).powi(2)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use epanet_rs::model::network::{JunctionData, Network, PipeData};
    use epanet_rs::model::options::HeadlossFormula;
    use epanet_rs::model::units::FlowUnits;

    fn sample_network() -> Network {
        let mut network = Network::new(FlowUnits::LPS, HeadlossFormula::DarcyWeisbach);
        network
            .add_junction(
                "A",
                &JunctionData {
                    elevation: 0.0,
                    coordinates: Some((0.0, 0.0)),
                    ..Default::default()
                },
            )
            .unwrap();
        network
            .add_junction(
                "B",
                &JunctionData {
                    elevation: 0.0,
                    coordinates: Some((100.0, 0.0)),
                    ..Default::default()
                },
            )
            .unwrap();
        network
            .add_pipe(
                "P1",
                &PipeData {
                    start_node: "A".into(),
                    end_node: "B".into(),
                    length: 100.0,
                    diameter: 12.0,
                    roughness: 100.0,
                    ..Default::default()
                },
            )
            .unwrap();
        network
    }

    #[test]
    fn finds_nearest_node() {
        let network = sample_network();
        let index = SpatialIndex::build(&network, 50.0);
        let found = index.nearest_node((2.0, 1.0), 10.0);
        assert_eq!(found, Some(network.node_map["A"]));
    }

    #[test]
    fn finds_nearest_link() {
        let network = sample_network();
        let index = SpatialIndex::build(&network, 50.0);
        let found = index.nearest_link((50.0, 1.0), 5.0);
        assert_eq!(found, Some(network.link_map["P1"]));
    }

    #[test]
    fn box_select_finds_nodes_within_bounds() {
        let network = sample_network();
        let index = SpatialIndex::build(&network, 50.0);
        let bounds = Bounds {
            min_x: -10.0,
            min_y: -10.0,
            max_x: 10.0,
            max_y: 10.0,
        };
        let found = index.nodes_in_bounds(bounds);
        assert_eq!(found, vec![network.node_map["A"]]);
    }
}
