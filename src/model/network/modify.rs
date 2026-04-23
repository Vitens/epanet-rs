use crate::error::InputError;
use crate::model::control::ControlCondition;
use crate::model::curve::{Curve, HeadCurve, ValveCurve};
use crate::model::junction::Junction;
use crate::model::pattern::Pattern;
use crate::model::link::{Link, LinkStatus, LinkType};
use crate::model::network::Network;
use crate::model::node::{Node, NodeType};
use crate::model::pipe::Pipe;
use crate::model::pump::Pump;
use crate::model::reservoir::Reservoir;
use crate::model::tank::Tank;
use crate::model::units::UnitConversion;
use crate::model::valve::{Valve, ValveType};

#[derive(Default)]
pub struct JunctionData {
  pub elevation: f64,
  pub basedemand: f64,
  pub emitter_coefficient: f64,
  pub pattern: Option<Box<str>>,
  pub coordinates: Option<(f64, f64)>,
}

#[derive(Default)]
pub struct JunctionUpdate {
  pub elevation: Option<f64>,
  pub basedemand: Option<f64>,
  pub emitter_coefficient: Option<f64>,
  /// Demand pattern id. `None` leaves the pattern untouched, `Some(None)`
  /// clears an existing pattern, `Some(Some(id))` sets the pattern.
  pub pattern: Option<Option<Box<str>>>,
  pub coordinates: Option<(f64, f64)>,
}

#[derive(Default)]
pub struct TankData {
  pub elevation: f64,
  pub initial_level: f64,
  pub min_level: f64,
  pub max_level: f64,
  pub diameter: f64,
  pub min_volume: f64,
  pub volume_curve_id: Option<Box<str>>,
  pub overflow: bool,
  pub coordinates: Option<(f64, f64)>,
}

#[derive(Default)]
pub struct TankUpdate {
  pub elevation: Option<f64>,
  pub coordinates: Option<(f64, f64)>,
  pub initial_level: Option<f64>,
  pub min_level: Option<f64>,
  pub max_level: Option<f64>,
  pub diameter: Option<f64>,
  pub min_volume: Option<f64>,
  pub overflow: Option<bool>,
}

/// Data for adding a new reservoir to the network.
#[derive(Default)]
pub struct ReservoirData {
  pub elevation: f64,
  pub head_pattern: Option<Box<str>>,
  pub coordinates: Option<(f64, f64)>,
}

/// Data for updating an existing reservoir.
#[derive(Default)]
pub struct ReservoirUpdate {
  pub elevation: Option<f64>,
  pub coordinates: Option<(f64, f64)>,
  /// Head pattern id. `None` leaves the pattern untouched, `Some(None)`
  /// clears an existing head pattern, `Some(Some(id))` sets the pattern.
  pub head_pattern: Option<Option<Box<str>>>,
}

/// Generic update for any node. Only elevation and coordinates can be changed
/// through this struct; type-specific properties must be updated via the
/// corresponding `update_tank`/`update_reservoir`/`update_junction` methods.
#[derive(Default)]
pub struct NodeUpdate {
  pub elevation: Option<f64>,
  pub coordinates: Option<(f64, f64)>,
}

/// Data for adding a pipe. All fields are in the user's unit system.
/// `minor_loss` is the dimensionless user-facing minor-loss coefficient (K).
/// The headloss formula is derived from the network's configured options.
pub struct PipeData {
  pub start_node: Box<str>,
  pub end_node: Box<str>,
  pub length: f64,
  pub diameter: f64,
  pub roughness: f64,
  pub minor_loss: f64,
  pub check_valve: bool,
  pub initial_status: LinkStatus,
  pub vertices: Option<Vec<(f64, f64)>>,
}

/// Data for updating a pipe.
/// All units are in the user's unit system.
#[derive(Default)]
pub struct PipeUpdate {
  pub length: Option<f64>,
  pub diameter: Option<f64>,
  pub roughness: Option<f64>,
  /// User-facing minor-loss coefficient (K). Internally re-normalized against
  /// the (possibly updated) diameter.
  pub minor_loss: Option<f64>,
  pub check_valve: Option<bool>,
}

/// Data for adding a pump. `head_curve_id` must reference an existing curve.
/// All units are in the user's unit system.
pub struct PumpData {
  pub start_node: Box<str>,
  pub end_node: Box<str>,
  pub speed: f64,
  pub head_curve_id: Option<Box<str>>,
  pub power: f64,
  pub initial_status: LinkStatus,
  pub vertices: Option<Vec<(f64, f64)>>,
}

/// Data for updating a pump.
/// All units are in the user's unit system.
#[derive(Default)]
pub struct PumpUpdate {
  pub speed: Option<f64>,
  pub power: Option<f64>,
  /// Head curve id. `None` leaves the curve untouched, `Some(None)` clears
  /// the existing head curve, `Some(Some(id))` sets the curve.
  pub head_curve_id: Option<Option<Box<str>>>,
}

/// Data for adding a valve. `curve_id` must reference an existing curve and is
/// only meaningful for GPV and PCV valves.
pub struct ValveData {
  pub start_node: Box<str>,
  pub end_node: Box<str>,
  pub diameter: f64,
  pub valve_type: ValveType,
  pub setting: f64,
  pub curve_id: Option<Box<str>>,
  pub minor_loss: f64,
  pub initial_status: LinkStatus,
  pub vertices: Option<Vec<(f64, f64)>>,
}

#[derive(Default)]
pub struct ValveUpdate {
  pub diameter: Option<f64>,
  pub setting: Option<f64>,
  pub minor_loss: Option<f64>,
  /// Valve curve id. `None` leaves the curve untouched, `Some(None)` clears
  /// the existing curve, `Some(Some(id))` sets the curve. GPV valves cannot
  /// have their curve cleared.
  pub curve_id: Option<Option<Box<str>>>,
  /// Change the valve's kind. Because the stored `setting` is interpreted
  /// differently per valve type, a new `setting` must be supplied whenever
  /// the type is changed. If the new type is `GPV`, a `curve_id` must also
  /// be supplied.
  pub valve_type: Option<ValveType>,
}

/// Generic update for any link. Changing `start_node`/`end_node` is a topology change, forcing a complete solver and state re-initialization.
/// changing `initial_status` or `vertices` is not.
#[derive(Default)]
pub struct LinkUpdate {
  pub start_node: Option<Box<str>>,
  pub end_node: Option<Box<str>>,
  pub vertices: Option<Vec<(f64, f64)>>,
  pub initial_status: Option<LinkStatus>,
}

/// Data for adding a new time pattern.
#[derive(Default)]
pub struct PatternData {
  pub multipliers: Vec<f64>,
}

/// Data for updating an existing pattern. Only provided fields are changed.
#[derive(Default)]
pub struct PatternUpdate {
  pub multipliers: Option<Vec<f64>>,
}

/// Data for adding a new curve. `x` and `y` must be non-empty and the same
/// length; `x` must be strictly increasing.
#[derive(Default)]
pub struct CurveData {
  pub x: Vec<f64>,
  pub y: Vec<f64>,
}

/// Data for updating an existing curve. Any provided axis replaces the
/// existing values; the resulting `x`/`y` pair must still satisfy the
/// non-empty / same-length / monotonic-x invariants.
#[derive(Default)]
pub struct CurveUpdate {
  pub x: Option<Vec<f64>>,
  pub y: Option<Vec<f64>>,
}

impl Network {
  /// Resets the change tracking flags by clearing the updated node and link sets.
  pub fn reset_changes(&mut self) {
    self.updated_nodes.clear();
    self.updated_links.clear();
  }

  /// Add a new junction to the network.
  pub fn add_junction(&mut self, id: &str, data: &JunctionData) -> Result<(), InputError> {

    // resolve the pattern index
    let pattern_index = if let Some(pattern) = &data.pattern {
      Some(self.pattern_map.get(pattern).ok_or(InputError::PatternNotFound { pattern_id: pattern.clone() })?)
    } else {
      None
    };

    let junction = Junction {
      basedemand: data.basedemand,
      emitter_coefficient: data.emitter_coefficient,
      pattern: data.pattern.clone(),
      pattern_index: pattern_index.copied(),
    };

    let mut node = Node {
      id: id.into(),
      elevation: data.elevation,
      node_type: NodeType::Junction(junction),
      coordinates: data.coordinates,
    };

    // convert the node to standard units
    node.convert_to_standard(&self.options);
    // add the node to the network
    self.add_node(node)?;

    Ok(())
  }

  /// Add a new tank to the network.
  pub fn add_tank(&mut self, id: &str, data: &TankData) -> Result<(), InputError> {
    // validate the volume curve id (if provided) exists
    if let Some(curve_id) = &data.volume_curve_id {
      if !self.curve_map.contains_key(curve_id) {
        return Err(InputError::CurveNotFound { curve_id: curve_id.clone() });
      }
    }

    let tank = Tank {
      elevation: data.elevation,
      initial_level: data.initial_level,
      min_level: data.min_level,
      max_level: data.max_level,
      diameter: data.diameter,
      min_volume: data.min_volume,
      volume_curve_id: data.volume_curve_id.clone(),
      overflow: data.overflow,
      volume_curve: None,
      links_to: Vec::new(),
      links_from: Vec::new(),
    };
    let mut node = Node {
      id: id.into(),
      elevation: data.elevation,
      node_type: NodeType::Tank(tank),
      coordinates: data.coordinates,
    };

    // convert the tank node to standard units
    node.convert_to_standard(&self.options);

    // add the node to the network
    self.add_node(node)?;

    Ok(())
  }

  /// Add a new reservoir to the network.
  pub fn add_reservoir(&mut self, id: &str, data: &ReservoirData) -> Result<(), InputError> {
    // resolve the head pattern index
    let head_pattern_index = if let Some(pattern) = &data.head_pattern {
      Some(*self.pattern_map.get(pattern).ok_or(InputError::PatternNotFound { pattern_id: pattern.clone() })?)
    } else {
      None
    };

    let reservoir = Reservoir {
      head_pattern: data.head_pattern.clone(),
      head_pattern_index,
    };
    let mut node = Node {
      id: id.into(),
      elevation: data.elevation,
      node_type: NodeType::Reservoir(reservoir),
      coordinates: data.coordinates,
    };

    // convert the reservoir node to standard units
    node.convert_to_standard(&self.options);
    self.add_node(node)?;

    Ok(())
  }

  pub fn update_junction(&mut self, id: &str, update: &JunctionUpdate) -> Result<(), InputError> {
    // lookup the node index
    let node_index = self.node_map.get(id).ok_or(InputError::NodeNotFound { node_id: id.into() })?;

    // resolve the new pattern index up front (if a pattern change is
    // requested) so we fail fast on an unknown pattern before mutating
    // anything. `Some(Some(id))` sets the pattern; `Some(None)` clears it.
    let pattern_change: Option<(Option<Box<str>>, Option<usize>)> = match &update.pattern {
      None => None, // do nothing
      Some(None) => Some((None, None)), // clear the pattern
      Some(Some(pattern_id)) => {
        let idx = *self.pattern_map.get(pattern_id)
          .ok_or_else(|| InputError::PatternNotFound { pattern_id: pattern_id.clone() })?;
        Some((Some(pattern_id.clone()), Some(idx))) // set the pattern and pattern index
      }
    };

    let node = &mut self.nodes[*node_index];
    // check if the node is a tank
    if !matches!(node.node_type, NodeType::Junction(_)) {
      return Err(InputError::NodeNotAJunction { node_id: id.into() });
    }

    // convert to user units
    node.convert_from_standard(&self.options);

    if let NodeType::Junction(junction) = &mut node.node_type {
      node.elevation = update.elevation.unwrap_or(node.elevation);

      if let Some(coordinates) = update.coordinates {
        node.coordinates = Some(coordinates);
      }

      junction.basedemand = update.basedemand.unwrap_or(junction.basedemand);

      junction.emitter_coefficient = update.emitter_coefficient.unwrap_or(junction.emitter_coefficient);

      if let Some((new_pattern, new_index)) = pattern_change {
        junction.pattern = new_pattern;
        junction.pattern_index = new_index;
      }
    }
    // convert back to standard units
    node.convert_to_standard(&self.options);

    // add the node index to the updated nodes set, forcing state update
    self.updated_nodes.insert(*node_index);
    // increment the properties version
    self.properties_version += 1;

    Ok(())
  }

  /// Update the elevation and/or coordinates of any node (junction, tank, or reservoir).
  /// For tanks, the tank's internal elevation is kept in sync with the node elevation.
  pub fn update_node(&mut self, id: &str, update: &NodeUpdate) -> Result<(), InputError> {
    let node_index = *self.node_map.get(id).ok_or(InputError::NodeNotFound { node_id: id.into() })?;

    let node = &mut self.nodes[node_index];

    if let Some(coordinates) = update.coordinates {
      node.coordinates = Some(coordinates);
    }

    if let Some(new_elevation) = update.elevation {
      node.convert_from_standard(&self.options);
      node.elevation = new_elevation;
      // if the node is a tank, update the tank's elevation as well
      if let NodeType::Tank(tank) = &mut node.node_type {
        tank.elevation = node.elevation;
      }
      node.convert_to_standard(&self.options);

      // elevation affects the initial head of fixed-head nodes, so mark it as updated
      self.updated_nodes.insert(node_index);
      self.properties_version += 1;
    }


    Ok(())
  }

  /// Update a tank's properties. All fields are optional; only provided fields are changed.
  pub fn update_tank(&mut self, id: &str, update: &TankUpdate) -> Result<(), InputError> {
    let node_index = *self.node_map.get(id).ok_or(InputError::NodeNotFound { node_id: id.into() })?;

    let node = &mut self.nodes[node_index];

    // check if the node is a tank
    if !matches!(node.node_type, NodeType::Tank(_)) {
      return Err(InputError::NodeNotATank { node_id: id.into() });
    }

    if let Some(coordinates) = update.coordinates {
      node.coordinates = Some(coordinates);
    }

    // convert to user units for in-place editing
    node.convert_from_standard(&self.options);

    if let Some(elevation) = update.elevation {
      node.elevation = elevation;
    }

    if let NodeType::Tank(tank) = &mut node.node_type {
      // keep tank.elevation in sync with node.elevation in user-unit space
      // so convert_to_standard below converts both consistently
      tank.elevation = node.elevation;
      tank.initial_level = update.initial_level.unwrap_or(tank.initial_level);
      tank.min_level = update.min_level.unwrap_or(tank.min_level);
      tank.max_level = update.max_level.unwrap_or(tank.max_level);
      tank.diameter = update.diameter.unwrap_or(tank.diameter);
      tank.min_volume = update.min_volume.unwrap_or(tank.min_volume);
      tank.overflow = update.overflow.unwrap_or(tank.overflow);
    }

    node.convert_to_standard(&self.options);

    self.updated_nodes.insert(node_index);
    self.properties_version += 1;

    Ok(())
  }

  /// Update a reservoir's properties. `head_pattern` follows the
  /// `Option<Option<_>>` convention: `None` leaves the pattern untouched,
  /// `Some(None)` clears it, `Some(Some(id))` sets it.
  pub fn update_reservoir(&mut self, id: &str, update: &ReservoirUpdate) -> Result<(), InputError> {
    let node_index = *self.node_map.get(id).ok_or(InputError::NodeNotFound { node_id: id.into() })?;

    // resolve the new pattern index up front (if a pattern change is
    // requested) so we fail fast on an unknown pattern.
    let pattern_change: Option<(Option<Box<str>>, Option<usize>)> = match &update.head_pattern {
      None => None, // do nothing
      Some(None) => Some((None, None)), // clear the pattern
      Some(Some(pattern_id)) => {
        let idx = *self.pattern_map.get(pattern_id)
          .ok_or_else(|| InputError::PatternNotFound { pattern_id: pattern_id.clone() })?;
        Some((Some(pattern_id.clone()), Some(idx))) // set the pattern and pattern index
      }
    };

    let node = &mut self.nodes[node_index];

    if !matches!(node.node_type, NodeType::Reservoir(_)) {
      return Err(InputError::NodeNotAReservoir { node_id: id.into() });
    }

    if let Some(coordinates) = update.coordinates {
      node.coordinates = Some(coordinates);
    }

    if let Some(elevation) = update.elevation {
      node.convert_from_standard(&self.options);
      node.elevation = elevation;
      node.convert_to_standard(&self.options);
    }

    if let NodeType::Reservoir(reservoir) = &mut node.node_type {
      if let Some((new_pattern, new_index)) = pattern_change {
        reservoir.head_pattern = new_pattern;
        reservoir.head_pattern_index = new_index;
      }
    }

    self.updated_nodes.insert(node_index);
    self.properties_version += 1;

    Ok(())
  }

  // --- Link add/update methods ---

  /// Add a new pipe to the network. Fails if either referenced node does not
  /// exist or a link with the same id is already present.
  pub fn add_pipe(&mut self, id: &str, data: &PipeData) -> Result<(), InputError> {
    let start_node = *self.node_map.get(&data.start_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.start_node.clone() })?;
    let end_node = *self.node_map.get(&data.end_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.end_node.clone() })?;

    let pipe = Pipe {
      diameter: data.diameter,
      length: data.length,
      roughness: data.roughness,
      minor_loss: data.minor_loss,
      check_valve: data.check_valve,
      headloss_formula: self.options.headloss_formula,
    };

    let mut link = Link {
      id: id.into(),
      link_type: LinkType::Pipe(pipe),
      start_node_id: data.start_node.clone(),
      end_node_id: data.end_node.clone(),
      initial_status: data.initial_status,
      vertices: data.vertices.clone(),
      start_node,
      end_node,
    };

    // convert the link to internal units (feet, etc.)
    link.convert_to_standard(&self.options);
    // normalize the user-facing minor-loss coefficient against the internal diameter
    if let LinkType::Pipe(pipe) = &mut link.link_type {
      pipe.minor_loss = 0.02517 * pipe.minor_loss / pipe.diameter.powi(4);
    }

    self.add_link(link)
  }

  /// Add a new pump to the network. Any `head_curve_id` is resolved eagerly.
  pub fn add_pump(&mut self, id: &str, data: &PumpData) -> Result<(), InputError> {
    let start_node = *self.node_map.get(&data.start_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.start_node.clone() })?;
    let end_node = *self.node_map.get(&data.end_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.end_node.clone() })?;

    let head_curve = self.resolve_head_curve(data.head_curve_id.as_deref())?;

    let pump = Pump {
      speed: data.speed,
      head_curve_id: data.head_curve_id.clone(),
      power: data.power,
      head_curve: None,
    };

    let mut link = Link {
      id: id.into(),
      link_type: LinkType::Pump(pump),
      start_node_id: data.start_node.clone(),
      end_node_id: data.end_node.clone(),
      initial_status: data.initial_status,
      vertices: data.vertices.clone(),
      start_node,
      end_node,
    };

    link.convert_to_standard(&self.options);
    // assign the (already standard-units) head curve after unit conversion so
    // it does not get double-converted by Pump::convert_to_standard
    if let LinkType::Pump(pump) = &mut link.link_type {
      pump.head_curve = head_curve;
    }

    self.add_link(link)
  }

  /// Add a new valve to the network. Any `curve_id` is resolved eagerly. For
  /// PSV/PRV valves the stored setting is offset by the attached node's
  /// elevation, matching the input file loader's behaviour.
  pub fn add_valve(&mut self, id: &str, data: &ValveData) -> Result<(), InputError> {
    let start_node = *self.node_map.get(&data.start_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.start_node.clone() })?;
    let end_node = *self.node_map.get(&data.end_node)
      .ok_or_else(|| InputError::NodeNotFound { node_id: data.end_node.clone() })?;

    // // GPV valves are driven by their curve; require it up front
    // if data.valve_type == ValveType::GPV && data.curve_id.is_none() {
    //   return Err(InputError::new("GPV valves require a curve id"));
    // }
    let valve_curve = self.resolve_valve_curve(data.curve_id.as_deref())?;

    let valve = Valve {
      diameter: data.diameter,
      setting: data.setting,
      curve_id: data.curve_id.clone(),
      valve_type: data.valve_type.clone(),
      minor_loss: data.minor_loss,
      valve_curve: None,
    };

    let mut link = Link {
      id: id.into(),
      link_type: LinkType::Valve(valve),
      start_node_id: data.start_node.clone(),
      end_node_id: data.end_node.clone(),
      initial_status: data.initial_status,
      vertices: data.vertices.clone(),
      start_node,
      end_node,
    };

    link.convert_to_standard(&self.options);

    // apply the node-elevation offset to PSV/PRV settings and attach the curve
    if let LinkType::Valve(valve) = &mut link.link_type {
      match valve.valve_type {
        ValveType::PSV => {
          valve.setting += self.nodes[start_node].elevation;
          self.contains_pressure_control_valve = true;
        }
        ValveType::PRV => {
          valve.setting += self.nodes[end_node].elevation;
          self.contains_pressure_control_valve = true;
        }
        _ => {}
      }
      valve.valve_curve = valve_curve;
    }

    self.add_link(link)
  }

  /// Update properties of an existing pipe. Only provided fields are changed.
  pub fn update_pipe(&mut self, id: &str, update: &PipeUpdate) -> Result<(), InputError> {
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })?;
    let link = &mut self.links[link_index];

    if !matches!(link.link_type, LinkType::Pipe(_)) {
      return Err(InputError::LinkNotAPipe { link_id: id.into() });
    }


    link.convert_from_standard(&self.options);

    if let LinkType::Pipe(pipe) = &mut link.link_type {
      pipe.length = update.length.unwrap_or(pipe.length);
      pipe.diameter = update.diameter.unwrap_or(pipe.diameter);
      pipe.roughness = update.roughness.unwrap_or(pipe.roughness);
      pipe.check_valve = update.check_valve.unwrap_or(pipe.check_valve);
    }

    link.convert_to_standard(&self.options);

    // re-normalize the minor-loss coefficient against the (possibly new) diameter
    if let LinkType::Pipe(pipe) = &mut link.link_type {
      if let Some(minor_loss) = update.minor_loss {
        // update the minor loss coefficient (this has to be done using standard units!)
        pipe.minor_loss = 0.02517 * minor_loss / pipe.diameter.powi(4);
      }
    }

    self.updated_links.insert(link_index);
    self.properties_version += 1;
    Ok(())
  }

  /// Update properties of an existing pump. `head_curve_id` follows the
  /// `Option<Option<_>>` convention: `None` leaves the curve untouched,
  /// `Some(None)` clears it (turning the pump into a constant-power or
  /// open-valve pump depending on `power`), and `Some(Some(id))` attaches
  /// the named head curve.
  pub fn update_pump(&mut self, id: &str, update: &PumpUpdate) -> Result<(), InputError> {
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })?;

    if !matches!(self.links[link_index].link_type, LinkType::Pump(_)) {
      return Err(InputError::LinkNotAPump { link_id: id.into() });
    }

    // resolve the new head curve (if a curve change is requested) before we
    // take a mutable borrow on the link.
    let curve_change: Option<(Option<Box<str>>, Option<HeadCurve>)> = match &update.head_curve_id {
      None => None, // do nothing
      Some(None) => Some((None, None)), // clear the curve
      Some(Some(curve_id)) => {
        let resolved = self.resolve_head_curve(Some(curve_id.as_ref()))?;
        Some((Some(curve_id.clone()), resolved)) // set the curve and curve index
      }
    };

    let link = &mut self.links[link_index];
    link.convert_from_standard(&self.options);

    if let LinkType::Pump(pump) = &mut link.link_type {
      if let Some(speed) = update.speed { pump.speed = speed; }
      if let Some(power) = update.power { pump.power = power; }
    }

    link.convert_to_standard(&self.options);

    if let LinkType::Pump(pump) = &mut link.link_type {
      if let Some((new_curve_id, new_head_curve)) = curve_change {
        pump.head_curve_id = new_curve_id;
        pump.head_curve = new_head_curve;
      }
    }

    self.updated_links.insert(link_index);
    self.properties_version += 1;
    Ok(())
  }

  /// Update properties of an existing valve. `curve_id` follows the
  /// `Option<Option<_>>` convention: `None` leaves the curve untouched,
  /// `Some(None)` clears it, `Some(Some(id))` sets it (GPV valves cannot
  /// have their curve cleared). The PSV/PRV elevation offset is
  /// automatically re-applied when the `setting` is changed. Changing
  /// `valve_type` requires a new `setting` (the stored value has no
  /// meaning under a different type); switching to `GPV` also requires a
  /// `curve_id` to be supplied.
  pub fn update_valve(&mut self, id: &str, update: &ValveUpdate) -> Result<(), InputError> {
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })?;

    if !matches!(self.links[link_index].link_type, LinkType::Valve(_)) {
      return Err(InputError::LinkNotAValve { link_id: id.into() });
    }

    // capture the current valve_type and resolve the target valve_type so we
    // can validate the update before mutating anything.
    let (old_valve_type, had_curve) = if let LinkType::Valve(valve) = &self.links[link_index].link_type {
      (valve.valve_type.clone(), valve.curve_id.is_some())
    } else {
      unreachable!()
    };
    let new_valve_type = update.valve_type.clone().unwrap_or_else(|| old_valve_type.clone());
    let type_changed = new_valve_type != old_valve_type;

    if type_changed && update.setting.is_none() {
      return Err(InputError::new(format!(
        "Changing valve_type on valve {} requires a new setting", id
      )));
    }

    // GPV valves always need a curve; the new rules for the tri-state
    // `curve_id` field are:
    //   * Some(None)         -> explicit clear is rejected for GPV
    //   * type_changed=true  -> caller must provide Some(Some(_))
    //   * type_changed=false -> keeping the existing curve (None) is fine
    if new_valve_type == ValveType::GPV {
      if matches!(update.curve_id, Some(None)) {
        return Err(InputError::new(format!(
          "GPV valve {} requires a curve id", id
        )));
      }
      if type_changed && !matches!(update.curve_id, Some(Some(_))) {
        return Err(InputError::new(format!(
          "GPV valve {} requires a curve id", id
        )));
      }
    }

    // resolve the new curve (if a curve change is requested) before we take
    // a mutable borrow on the link.
    let curve_change: Option<(Option<Box<str>>, Option<ValveCurve>)> = match &update.curve_id {
      None => None,
      Some(None) => Some((None, None)),
      Some(Some(curve_id)) => {
        let resolved = self.resolve_valve_curve(Some(curve_id.as_ref()))?;
        Some((Some(curve_id.clone()), resolved))
      }
    };

    // capture start/end node elevations (already in standard units) for PSV/PRV
    let start_elevation = self.nodes[self.links[link_index].start_node].elevation;
    let end_elevation = self.nodes[self.links[link_index].end_node].elevation;

    let link = &mut self.links[link_index];

    let stored_setting = if let LinkType::Valve(valve) = &link.link_type {
      valve.setting
    } else {
      unreachable!()
    };

    // strip the old PSV/PRV elevation offset before converting to user units.
    // The offset for the (possibly new) type is re-applied after the
    // to-standard step below, matching add_valve's semantics.
    if let LinkType::Valve(valve) = &mut link.link_type {
      match old_valve_type {
        ValveType::PSV => { valve.setting = stored_setting - start_elevation; }
        ValveType::PRV => { valve.setting = stored_setting - end_elevation; }
        _ => {}
      }
    }

    // convert_from_standard uses valve.valve_type, which is still the OLD
    // type here so the existing `setting` is interpreted correctly.
    link.convert_from_standard(&self.options);

    if let LinkType::Valve(valve) = &mut link.link_type {
      if let Some(diameter) = update.diameter { valve.diameter = diameter; }
      if let Some(setting) = update.setting { valve.setting = setting; }
      if let Some(minor_loss) = update.minor_loss { valve.minor_loss = minor_loss; }
      // switch to the NEW type before converting back to standard so that
      // the new `setting` is interpreted in new-type semantics.
      valve.valve_type = new_valve_type.clone();
    }

    link.convert_to_standard(&self.options);

    if let LinkType::Valve(valve) = &mut link.link_type {
      match new_valve_type {
        ValveType::PSV => { valve.setting += start_elevation; }
        ValveType::PRV => { valve.setting += end_elevation; }
        _ => {}
      }
      if let Some((new_curve_id, new_valve_curve)) = curve_change {
        valve.curve_id = new_curve_id;
        valve.valve_curve = new_valve_curve;
      }
    }

    // Recompute the PSV/PRV flag when the type changed (we may have added a
    // pressure-control valve or removed the last one). Otherwise just make
    // sure the flag is set for existing PSV/PRV valves.
    if type_changed {
      self.contains_pressure_control_valve = self.links.iter().any(|l| {
        matches!(&l.link_type, LinkType::Valve(v)
          if matches!(v.valve_type, ValveType::PSV | ValveType::PRV))
      });
    } else if matches!(new_valve_type, ValveType::PSV | ValveType::PRV) {
      self.contains_pressure_control_valve = true;
    }

    self.updated_links.insert(link_index);
    self.properties_version += 1;
    Ok(())
  }

  /// Update the topology/display properties of any link: endpoints, vertices
  /// and initial status. Type-specific properties must be edited via
  /// [`update_pipe`], [`update_pump`] or [`update_valve`].
  pub fn update_link(&mut self, id: &str, update: &LinkUpdate) -> Result<(), InputError> {
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })?;

    // resolve new endpoints up front so we fail fast on missing nodes
    let new_start = if let Some(start_id) = &update.start_node {
      Some((start_id.clone(), *self.node_map.get(start_id)
        .ok_or_else(|| InputError::NodeNotFound { node_id: start_id.clone() })?))
    } else {
      None
    };
    let new_end = if let Some(end_id) = &update.end_node {
      Some((end_id.clone(), *self.node_map.get(end_id)
        .ok_or_else(|| InputError::NodeNotFound { node_id: end_id.clone() })?))
    } else {
      None
    };

    let link = &self.links[link_index];
    let old_start = link.start_node;
    let old_end = link.end_node;

    // remember the valve_type so we can re-offset PSV/PRV settings if the
    // attached node (and thus its elevation) is swapped out
    let valve_type = if let LinkType::Valve(valve) = &link.link_type {
      Some(valve.valve_type.clone())
    } else {
      None
    };

    let mut topology_changed = false;

    if let Some((start_id, start_idx)) = new_start {
      let old_elev = self.nodes[old_start].elevation;
      let new_elev = self.nodes[start_idx].elevation;

      // detach from old tank's links_from and attach to new tank's links_from
      remove_tank_link(&mut self.nodes[old_start], link_index, /*outgoing=*/true);
      attach_tank_link(&mut self.nodes[start_idx], link_index, /*outgoing=*/true);

      let link = &mut self.links[link_index];
      link.start_node = start_idx;
      link.start_node_id = start_id;

      // keep PSV setting offset in sync with its start-node elevation
      if valve_type.as_ref() == Some(&ValveType::PSV) {
        if let LinkType::Valve(valve) = &mut link.link_type {
          valve.setting += new_elev - old_elev;
        }
      }
      topology_changed = true;
    }

    if let Some((end_id, end_idx)) = new_end {
      let old_elev = self.nodes[old_end].elevation;
      let new_elev = self.nodes[end_idx].elevation;

      remove_tank_link(&mut self.nodes[old_end], link_index, /*outgoing=*/false);
      attach_tank_link(&mut self.nodes[end_idx], link_index, /*outgoing=*/false);

      let link = &mut self.links[link_index];
      link.end_node = end_idx;
      link.end_node_id = end_id;

      // keep PRV setting offset in sync with its end-node elevation
      if valve_type.as_ref() == Some(&ValveType::PRV) {
        if let LinkType::Valve(valve) = &mut link.link_type {
          valve.setting += new_elev - old_elev;
        }
      }
      topology_changed = true;
    }

    let link = &mut self.links[link_index];
    if let Some(vertices) = &update.vertices {
      link.vertices = Some(vertices.clone());
    }

    let mut properties_changed = false;
    if let Some(status) = update.initial_status {
      if link.initial_status != status {
        link.initial_status = status;
        properties_changed = true;
      }
    }

    if topology_changed {
      self.topology_version += 1;
    }
    if properties_changed {
      self.updated_links.insert(link_index);
      self.properties_version += 1;
    }

    Ok(())
  }

  // --- removal methods ---

  /// Remove a link from the network. Fails if any control references the
  /// link. Uses `swap_remove` internally, so the last link in `links` takes
  /// the removed link's slot; `link_map`, tank link lists and `updated_links`
  /// are all kept in sync.
  /// `unconditional` is a boolean flag that indicates if the link should be removed if it is referenced by a control. If set to true, remove the link as well as any controls that reference it.
  pub fn remove_link(&mut self, id: &str, unconditional: bool) -> Result<(), InputError> {
    let link_index = *self.link_map.get(id)
      .ok_or_else(|| InputError::LinkNotFound { link_id: id.into() })?;

    if !unconditional && self.controls.iter().any(|c| c.link_id.as_ref() == id) {
      return Err(InputError::new(format!(
        "Cannot remove link {}: referenced by a control", id
      )));
    }

    let (start_node, end_node, was_pressure_control) = {
      let link = &self.links[link_index];
      let was_pc = matches!(&link.link_type, LinkType::Valve(v)
        if matches!(v.valve_type, ValveType::PSV | ValveType::PRV));
      (link.start_node, link.end_node, was_pc)
    };

    // detach from tank link lists at both endpoints
    remove_tank_link(&mut self.nodes[start_node], link_index, /*outgoing=*/true);
    remove_tank_link(&mut self.nodes[end_node], link_index, /*outgoing=*/false);

    self.link_map.remove(id);
    self.updated_links.remove(&link_index);

    let last_index = self.links.len() - 1;
    self.links.swap_remove(link_index);

    // if swap_remove moved a different link into `link_index`, reindex it
    if link_index != last_index {
      let (moved_id, moved_start, moved_end) = {
        let moved = &self.links[link_index];
        (moved.id.clone(), moved.start_node, moved.end_node)
      };
      self.link_map.insert(moved_id, link_index);
      retarget_tank_link(&mut self.nodes[moved_start], last_index, link_index, /*outgoing=*/true);
      retarget_tank_link(&mut self.nodes[moved_end], last_index, link_index, /*outgoing=*/false);
      if self.updated_links.remove(&last_index) {
        self.updated_links.insert(link_index);
      }
    }

    if was_pressure_control {
      self.contains_pressure_control_valve = self.links.iter().any(|l| {
        matches!(&l.link_type, LinkType::Valve(v)
          if matches!(v.valve_type, ValveType::PSV | ValveType::PRV))
      });
    }
    // remove any controls that reference this link
    if unconditional {
      self.controls.retain(|c| c.link_id.as_ref() != id);
    }

    self.topology_version += 1;
    Ok(())
  }

  /// Remove a node from the network. Any links attached to the node are
  /// cascade-removed first via `remove_link` (so their controls must not
  /// block removal). Fails if a control still references the node itself.
  /// Uses `swap_remove` internally, so the last node in `nodes` takes the
  /// removed node's slot; `node_map`, link endpoint indices, control
  /// node/tank indices and `updated_nodes` are all kept in sync.
  /// `unconditional` is a boolean flag that indicates if the node should be removed if it is referenced by a control. If set to true, remove the node as well as any controls that reference it.
  pub fn remove_node(&mut self, id: &str, unconditional: bool) -> Result<(), InputError> {
    let node_index = *self.node_map.get(id)
      .ok_or_else(|| InputError::NodeNotFound { node_id: id.into() })?;

    let referenced_by_control = self.controls.iter().any(|c| match &c.condition {
      ControlCondition::HighPressure { node_index: ni, .. }
      | ControlCondition::LowPressure { node_index: ni, .. } => *ni == node_index,
      ControlCondition::HighLevel { tank_index, .. }
      | ControlCondition::LowLevel { tank_index, .. } => *tank_index == node_index,
      _ => false,
    });
    if !unconditional && referenced_by_control {
      return Err(InputError::new(format!(
        "Cannot remove node {}: referenced by a control", id
      )));
    }

    // cascade-remove every link attached to this node (at either endpoint).
    // ids are collected up front since `remove_link` uses `swap_remove` and
    // therefore invalidates positional link indices as it runs.
    let connected_link_ids: Vec<Box<str>> = self.links.iter()
      .filter(|l| l.start_node == node_index || l.end_node == node_index)
      .map(|l| l.id.clone())
      .collect();

    // check if any of the connected links are referenced by a control
    if !unconditional && self.controls.iter().any(|c| connected_link_ids.contains(&c.link_id)) {
      return Err(InputError::new(format!(
        "Cannot remove node {}: referenced by a control", id
      )));
    }

    for link_id in &connected_link_ids {
      self.remove_link(link_id, unconditional)?;
    }

    // remove any controls that reference this node (by node or tank index).
    // must happen before `swap_remove` so the match is against the current
    // `node_index` rather than the post-swap slot.
    if unconditional {
      self.controls.retain(|c| match &c.condition {
        ControlCondition::HighPressure { node_index: ni, .. }
        | ControlCondition::LowPressure { node_index: ni, .. } => *ni != node_index,
        ControlCondition::HighLevel { tank_index, .. }
        | ControlCondition::LowLevel { tank_index, .. } => *tank_index != node_index,
        _ => true,
      });
    }

    self.node_map.remove(id);
    self.updated_nodes.remove(&node_index);

    let last_index = self.nodes.len() - 1;
    self.nodes.swap_remove(node_index);

    if node_index != last_index {
      let moved_id = self.nodes[node_index].id.clone();
      self.node_map.insert(moved_id, node_index);

      for link in self.links.iter_mut() {
        if link.start_node == last_index { link.start_node = node_index; }
        if link.end_node == last_index { link.end_node = node_index; }
      }

      for control in self.controls.iter_mut() {
        match &mut control.condition {
          ControlCondition::HighPressure { node_index: ni, .. }
          | ControlCondition::LowPressure { node_index: ni, .. } => {
            if *ni == last_index { *ni = node_index; }
          }
          ControlCondition::HighLevel { tank_index, .. }
          | ControlCondition::LowLevel { tank_index, .. } => {
            if *tank_index == last_index { *tank_index = node_index; }
          }
          _ => {}
        }
      }

      if self.updated_nodes.remove(&last_index) {
        self.updated_nodes.insert(node_index);
      }
    }

    self.topology_version += 1;
    Ok(())
  }

  // --- pattern methods ---

  /// Add a new time pattern to the network. Fails if a pattern with the same
  /// id already exists.
  pub fn add_pattern(&mut self, id: &str, data: &PatternData) -> Result<(), InputError> {
    if self.pattern_map.contains_key(id) {
      return Err(InputError::new(format!("Pattern {} already exists", id)));
    }
    let pattern = Pattern {
      id: id.into(),
      multipliers: data.multipliers.clone(),
    };
    self.pattern_map.insert(pattern.id.clone(), self.patterns.len());
    self.patterns.push(pattern);
    Ok(())
  }

  /// Update an existing pattern's multipliers in place. Only provided fields
  /// are changed.
  pub fn update_pattern(&mut self, id: &str, update: &PatternUpdate) -> Result<(), InputError> {
    let pattern_index = *self.pattern_map.get(id)
      .ok_or_else(|| InputError::PatternNotFound { pattern_id: id.into() })?;

    if let Some(multipliers) = &update.multipliers {
      self.patterns[pattern_index].multipliers = multipliers.clone();
      // pattern coefficients are re-read every timestep, so no per-node
      // invalidation is required; bump properties_version so callers can
      // still detect a change.
      self.properties_version += 1;
    }
    Ok(())
  }

  /// Remove a pattern from the network. Fails if any junction or reservoir
  /// still references the pattern. Uses `swap_remove` internally, so the last
  /// pattern in `patterns` takes the removed pattern's slot; `pattern_map`
  /// and every cached `pattern_index`/`head_pattern_index` is kept in sync.
  pub fn remove_pattern(&mut self, id: &str) -> Result<(), InputError> {
    let pattern_index = *self.pattern_map.get(id)
      .ok_or_else(|| InputError::PatternNotFound { pattern_id: id.into() })?;

    let referenced = self.nodes.iter().any(|n| match &n.node_type {
      NodeType::Junction(j) => j.pattern.as_deref() == Some(id),
      NodeType::Reservoir(r) => r.head_pattern.as_deref() == Some(id),
      _ => false,
    });
    if referenced {
      return Err(InputError::new(format!(
        "Cannot remove pattern {}: referenced by one or more nodes", id
      )));
    }

    self.pattern_map.remove(id);

    let last_index = self.patterns.len() - 1;
    self.patterns.swap_remove(pattern_index);

    // if swap_remove moved a different pattern into this slot, rebind it in
    // the map and re-resolve every node's cached pattern_index so stale
    // `last_index` hits become `pattern_index`.
    if pattern_index != last_index {
      let moved_id = self.patterns[pattern_index].id.clone();
      self.pattern_map.insert(moved_id, pattern_index);
      self.resolve_pattern_indices();
    }
    Ok(())
  }

  // --- curve methods ---

  /// Add a new curve to the network. Fails if the id already exists or the
  /// x/y arrays are empty, of unequal length, or x is not strictly
  /// increasing.
  pub fn add_curve(&mut self, id: &str, data: &CurveData) -> Result<(), InputError> {
    if self.curve_map.contains_key(id) {
      return Err(InputError::new(format!("Curve {} already exists", id)));
    }
    validate_curve_axes(&data.x, &data.y)?;
    let curve = Curve {
      id: id.into(),
      x: data.x.clone(),
      y: data.y.clone(),
    };
    self.curve_map.insert(curve.id.clone(), self.curves.len());
    self.curves.push(curve);
    Ok(())
  }

  /// Update an existing curve's axes. Any axis left `None` is preserved.
  /// Every pump/valve whose cached `HeadCurve`/`ValveCurve` was derived from
  /// this curve is rebuilt, and those links are marked as updated.
  pub fn update_curve(&mut self, id: &str, update: &CurveUpdate) -> Result<(), InputError> {
    let curve_index = *self.curve_map.get(id)
      .ok_or_else(|| InputError::CurveNotFound { curve_id: id.into() })?;

    // decide the new x/y up front so we can validate before mutating
    let new_x = update.x.clone().unwrap_or_else(|| self.curves[curve_index].x.clone());
    let new_y = update.y.clone().unwrap_or_else(|| self.curves[curve_index].y.clone());
    validate_curve_axes(&new_x, &new_y)?;

    self.curves[curve_index].x = new_x;
    self.curves[curve_index].y = new_y;

    // rebuild every derived HeadCurve/ValveCurve that references this curve
    // by id, and mark the owning link as updated so the solver picks up the
    // new coefficients.
    let curve = &self.curves[curve_index];
    for (i, link) in self.links.iter_mut().enumerate() {
      match &mut link.link_type {
        LinkType::Pump(pump) if pump.head_curve_id.as_deref() == Some(id) => {
          pump.head_curve = Some(HeadCurve::new(curve, &self.options.flow_units, &self.options.unit_system)?);
          self.updated_links.insert(i);
        }
        LinkType::Valve(valve) if valve.curve_id.as_deref() == Some(id) => {
          valve.valve_curve = Some(ValveCurve::new(curve, &self.options.flow_units, &self.options.unit_system)?);
          self.updated_links.insert(i);
        }
        _ => {}
      }
    }
    self.properties_version += 1;
    Ok(())
  }

  /// Remove a curve from the network. Fails if any pump, valve or tank still
  /// references the curve. Uses `swap_remove` internally; `curve_map` is
  /// kept in sync. Derived `HeadCurve`/`ValveCurve` caches are not tracked
  /// by index and therefore need no remapping.
  pub fn remove_curve(&mut self, id: &str) -> Result<(), InputError> {
    let curve_index = *self.curve_map.get(id)
      .ok_or_else(|| InputError::CurveNotFound { curve_id: id.into() })?;

    let referenced_by_link = self.links.iter().any(|l| match &l.link_type {
      LinkType::Pump(p) => p.head_curve_id.as_deref() == Some(id),
      LinkType::Valve(v) => v.curve_id.as_deref() == Some(id),
      _ => false,
    });
    let referenced_by_tank = self.nodes.iter().any(|n| match &n.node_type {
      NodeType::Tank(t) => t.volume_curve_id.as_deref() == Some(id),
      _ => false,
    });
    if referenced_by_link || referenced_by_tank {
      return Err(InputError::new(format!(
        "Cannot remove curve {}: referenced by one or more links or tanks", id
      )));
    }

    self.curve_map.remove(id);

    let last_index = self.curves.len() - 1;
    self.curves.swap_remove(curve_index);

    if curve_index != last_index {
      let moved_id = self.curves[curve_index].id.clone();
      self.curve_map.insert(moved_id, curve_index);
    }
    Ok(())
  }

  // --- private helpers ---

  fn resolve_head_curve(&self, curve_id: Option<&str>) -> Result<Option<HeadCurve>, InputError> {
    let Some(curve_id) = curve_id else { return Ok(None) };
    let curve_index = *self.curve_map.get(curve_id)
      .ok_or_else(|| InputError::CurveNotFound { curve_id: curve_id.into() })?;
    let curve = &self.curves[curve_index];
    Ok(Some(HeadCurve::new(curve, &self.options.flow_units, &self.options.unit_system)?))
  }

  fn resolve_valve_curve(&self, curve_id: Option<&str>) -> Result<Option<ValveCurve>, InputError> {
    let Some(curve_id) = curve_id else { return Ok(None) };
    let curve_index = *self.curve_map.get(curve_id)
      .ok_or_else(|| InputError::CurveNotFound { curve_id: curve_id.into() })?;
    let curve = &self.curves[curve_index];
    Ok(Some(ValveCurve::new(curve, &self.options.flow_units, &self.options.unit_system)?))
  }
}

/// Validate a raw curve's x/y arrays: both must be non-empty, of equal
/// length, and `x` must be strictly increasing. Unit-specific validation
/// (e.g. pump curve monotonicity) happens later when the derived
/// `HeadCurve`/`ValveCurve` is built.
fn validate_curve_axes(x: &[f64], y: &[f64]) -> Result<(), InputError> {
  if x.is_empty() || y.is_empty() {
    return Err(InputError::new("Curve x/y arrays must be non-empty"));
  }
  if x.len() != y.len() {
    return Err(InputError::new(format!(
      "Curve x/y arrays must have the same length (got {} and {})",
      x.len(), y.len(),
    )));
  }
  for i in 1..x.len() {
    if x[i] <= x[i - 1] {
      return Err(InputError::new("Curve x values must be strictly increasing"));
    }
  }
  Ok(())
}

/// Remove `link_index` from the tank's outgoing (`links_from`) or incoming
/// (`links_to`) list when the node happens to be a tank. No-op for other node
/// types.
fn remove_tank_link(node: &mut Node, link_index: usize, outgoing: bool) {
  if let NodeType::Tank(tank) = &mut node.node_type {
    let list = if outgoing { &mut tank.links_from } else { &mut tank.links_to };
    if let Some(pos) = list.iter().position(|&i| i == link_index) {
      list.swap_remove(pos);
    }
  }
}

/// Append `link_index` to the tank's outgoing or incoming list when the node
/// happens to be a tank. No-op for other node types.
fn attach_tank_link(node: &mut Node, link_index: usize, outgoing: bool) {
  if let NodeType::Tank(tank) = &mut node.node_type {
    let list = if outgoing { &mut tank.links_from } else { &mut tank.links_to };
    list.push(link_index);
  }
}

/// Replace `old_index` with `new_index` in the tank's outgoing (`links_from`)
/// or incoming (`links_to`) list. Used after `swap_remove` on `links` to
/// point tank link lists at the moved link's new slot. No-op for non-tank
/// nodes.
fn retarget_tank_link(node: &mut Node, old_index: usize, new_index: usize, outgoing: bool) {
  if let NodeType::Tank(tank) = &mut node.node_type {
    let list = if outgoing { &mut tank.links_from } else { &mut tank.links_to };
    for i in list.iter_mut() {
      if *i == old_index { *i = new_index; }
    }
  }
}


#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::MperFT;
  use crate::model::options::HeadlossFormula;
  use crate::model::reservoir::Reservoir;
  use crate::model::tank::Tank;
  use crate::model::units::FlowUnits;

  // --- helpers for tank/reservoir tests (no public add_tank/add_reservoir yet) ---

  fn test_tank_node(id: &str, elevation: f64) -> Node {
    Node {
      id: id.into(),
      elevation,
      node_type: NodeType::Tank(Tank {
        elevation,
        initial_level: 10.0,
        min_level: 0.0,
        max_level: 20.0,
        diameter: 50.0,
        min_volume: 0.0,
        volume_curve_id: None,
        overflow: false,
        volume_curve: None,
        links_to: Vec::new(),
        links_from: Vec::new(),
      }),
      coordinates: None,
    }
  }

  fn test_reservoir_node(id: &str, elevation: f64) -> Node {
    Node {
      id: id.into(),
      elevation,
      node_type: NodeType::Reservoir(Reservoir {
        head_pattern: None,
        head_pattern_index: None,
      }),
      coordinates: None,
    }
  }

  // --- Junction tests ---

  #[test]
  fn test_add_and_update_junction() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();

    assert_eq!(network.nodes.len(), 1);
    assert_eq!(network.nodes[0].id, "J1".into());
    assert_eq!(network.nodes[0].elevation, 100.0);
    assert_eq!(network.topology_version, 1);
    assert_eq!(network.properties_version, 0);

    let update = JunctionUpdate {
      elevation: Some(200.0),
      basedemand: Some(20.0),
      emitter_coefficient: Some(0.5),
      pattern: None,
      coordinates: None,
    };
    network.update_junction("J1", &update).unwrap();
    assert_eq!(network.nodes[0].elevation, 200.0);
    let NodeType::Junction(junction) = &network.nodes[0].node_type else {
      panic!("Expected Junction node type");
    };
    assert_eq!(junction.basedemand, 20.0);
    assert!((junction.emitter_coefficient - 9.231479).abs() < 1e-6);
    assert_eq!(network.properties_version, 1);
  }

  #[test]
  fn test_update_junction_none_pattern_leaves_existing_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0] });
    network.pattern_map.insert("P1".into(), 0);
    network.add_junction("J1", &JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: Some("P1".into()),
      coordinates: None,
    }).unwrap();

    network.update_junction("J1", &JunctionUpdate {
      basedemand: Some(20.0),
      pattern: None,
      ..Default::default()
    }).unwrap();

    let NodeType::Junction(junction) = &network.nodes[0].node_type else {
      panic!("Expected Junction node type");
    };
    assert_eq!(junction.pattern.as_deref(), Some("P1"));
    assert_eq!(junction.pattern_index, Some(0));
    assert_eq!(junction.basedemand, 20.0);
  }

  #[test]
  fn test_update_junction_clears_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0] });
    network.pattern_map.insert("P1".into(), 0);
    network.add_junction("J1", &JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: Some("P1".into()),
      coordinates: None,
    }).unwrap();

    network.update_junction("J1", &JunctionUpdate {
      pattern: Some(None),
      ..Default::default()
    }).unwrap();

    let NodeType::Junction(junction) = &network.nodes[0].node_type else {
      panic!("Expected Junction node type");
    };
    assert!(junction.pattern.is_none());
    assert!(junction.pattern_index.is_none());
  }

  #[test]
  fn test_update_junction_sets_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0] });
    network.pattern_map.insert("P1".into(), 0);
    network.add_junction("J1", &JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    }).unwrap();

    network.update_junction("J1", &JunctionUpdate {
      pattern: Some(Some("P1".into())),
      ..Default::default()
    }).unwrap();

    let NodeType::Junction(junction) = &network.nodes[0].node_type else {
      panic!("Expected Junction node type");
    };
    assert_eq!(junction.pattern.as_deref(), Some("P1"));
    assert_eq!(junction.pattern_index, Some(0));
  }

  #[test]
  fn test_update_junction_pattern_not_found() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    }).unwrap();

    let err = network.update_junction("J1", &JunctionUpdate {
      pattern: Some(Some("missing".into())),
      ..Default::default()
    }).unwrap_err();
    assert!(matches!(err, InputError::PatternNotFound { .. }));
  }

  #[test]
  fn test_add_junction_si() {
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 10.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();
    assert_eq!(network.nodes[0].elevation, 100.0 / MperFT);
  }

  // --- add_tank tests ---

  fn test_tank_data(elevation: f64) -> TankData {
    TankData {
      elevation,
      initial_level: 10.0,
      min_level: 0.0,
      max_level: 20.0,
      diameter: 50.0,
      min_volume: 0.0,
      volume_curve_id: None,
      overflow: false,
      coordinates: None,
    }
  }

  #[test]
  fn test_add_tank() {
    let mut network = Network::default();
    let data = TankData {
      coordinates: Some((1.0, 2.0)),
      ..test_tank_data(100.0)
    };
    network.add_tank("T1", &data).unwrap();

    assert_eq!(network.nodes.len(), 1);
    assert_eq!(network.nodes[0].id, "T1".into());
    assert_eq!(network.nodes[0].elevation, 100.0);
    assert_eq!(network.nodes[0].coordinates, Some((1.0, 2.0)));
    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.elevation, 100.0);
    assert_eq!(tank.initial_level, 10.0);
    assert_eq!(tank.max_level, 20.0);
    assert_eq!(tank.diameter, 50.0);
    assert_eq!(network.topology_version, 1);
    assert_eq!(network.node_map.get("T1").copied(), Some(0));
  }

  #[test]
  fn test_add_tank_si_units() {
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    network.add_tank("T1", &test_tank_data(100.0)).unwrap();

    // user-provided meters should be stored as feet internally; node and tank
    // elevations should stay in sync
    assert!((network.nodes[0].elevation - 100.0 / MperFT).abs() < 1e-9);
    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert!((tank.elevation - 100.0 / MperFT).abs() < 1e-9);
    assert!((tank.initial_level - 10.0 / MperFT).abs() < 1e-9);
    assert!((tank.diameter - 50.0 / MperFT).abs() < 1e-9);
  }

  #[test]
  fn test_add_tank_duplicate_id() {
    let mut network = Network::default();
    network.add_tank("T1", &test_tank_data(100.0)).unwrap();
    let err = network.add_tank("T1", &test_tank_data(50.0)).unwrap_err();
    assert!(matches!(err, InputError::NodeExists { .. }));
  }

  #[test]
  fn test_add_tank_unknown_volume_curve() {
    let mut network = Network::default();
    let data = TankData {
      volume_curve_id: Some("missing".into()),
      ..test_tank_data(100.0)
    };
    let err = network.add_tank("T1", &data).unwrap_err();
    assert!(matches!(err, InputError::CurveNotFound { .. }));
    assert_eq!(network.nodes.len(), 0);
    assert_eq!(network.topology_version, 0);
  }

  // --- add_reservoir tests ---

  #[test]
  fn test_add_reservoir() {
    let mut network = Network::default();
    let data = ReservoirData {
      elevation: 200.0,
      head_pattern: None,
      coordinates: Some((5.0, 6.0)),
    };
    network.add_reservoir("R1", &data).unwrap();

    assert_eq!(network.nodes.len(), 1);
    assert_eq!(network.nodes[0].id, "R1".into());
    assert_eq!(network.nodes[0].elevation, 200.0);
    assert_eq!(network.nodes[0].coordinates, Some((5.0, 6.0)));
    let NodeType::Reservoir(reservoir) = &network.nodes[0].node_type else {
      panic!("Expected Reservoir node type");
    };
    assert!(reservoir.head_pattern.is_none());
    assert!(reservoir.head_pattern_index.is_none());
    assert_eq!(network.topology_version, 1);
  }

  #[test]
  fn test_add_reservoir_with_head_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0, 2.0] });
    network.pattern_map.insert("P1".into(), 0);

    let data = ReservoirData {
      elevation: 200.0,
      head_pattern: Some("P1".into()),
      coordinates: None,
    };
    network.add_reservoir("R1", &data).unwrap();

    let NodeType::Reservoir(reservoir) = &network.nodes[0].node_type else {
      panic!("Expected Reservoir node type");
    };
    assert_eq!(reservoir.head_pattern.as_deref(), Some("P1"));
    assert_eq!(reservoir.head_pattern_index, Some(0));
  }

  #[test]
  fn test_add_reservoir_si_units() {
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    let data = ReservoirData {
      elevation: 100.0,
      head_pattern: None,
      coordinates: None,
    };
    network.add_reservoir("R1", &data).unwrap();
    assert!((network.nodes[0].elevation - 100.0 / MperFT).abs() < 1e-9);
  }

  #[test]
  fn test_add_reservoir_unknown_pattern() {
    let mut network = Network::default();
    let data = ReservoirData {
      elevation: 200.0,
      head_pattern: Some("missing".into()),
      coordinates: None,
    };
    let err = network.add_reservoir("R1", &data).unwrap_err();
    assert!(matches!(err, InputError::PatternNotFound { .. }));
    assert_eq!(network.nodes.len(), 0);
    assert_eq!(network.topology_version, 0);
  }

  #[test]
  fn test_add_reservoir_duplicate_id() {
    let mut network = Network::default();
    network.add_reservoir("R1", &ReservoirData {
      elevation: 100.0,
      head_pattern: None,
      coordinates: None,
    }).unwrap();
    let err = network.add_reservoir("R1", &ReservoirData {
      elevation: 50.0,
      head_pattern: None,
      coordinates: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeExists { .. }));
  }

  // --- NodeUpdate tests ---

  #[test]
  fn test_update_node_coordinates_only() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();
    let props_before = network.properties_version;

    network.update_node("J1", &NodeUpdate {
      elevation: None,
      coordinates: Some((1.0, 2.0)),
    }).unwrap();

    assert_eq!(network.nodes[0].coordinates, Some((1.0, 2.0)));
    assert_eq!(network.nodes[0].elevation, 100.0);
    // coordinates are not a solver-relevant property; properties_version should not change
    assert_eq!(network.properties_version, props_before);
    assert!(network.updated_nodes.is_empty());
  }

  #[test]
  fn test_update_node_elevation() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 100.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();

    network.update_node("J1", &NodeUpdate {
      elevation: Some(250.0),
      coordinates: None,
    }).unwrap();

    assert_eq!(network.nodes[0].elevation, 250.0);
    assert_eq!(network.properties_version, 1);
    assert!(network.updated_nodes.contains(&0));
  }

  #[test]
  fn test_update_node_elevation_keeps_tank_elevation_in_sync() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();

    network.update_node("T1", &NodeUpdate {
      elevation: Some(150.0),
      coordinates: None,
    }).unwrap();

    assert_eq!(network.nodes[0].elevation, 150.0);
    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.elevation, 150.0);
  }

  #[test]
  fn test_update_node_not_found() {
    let mut network = Network::default();
    let err = network.update_node("missing", &NodeUpdate {
      elevation: Some(1.0),
      coordinates: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
  }

  // --- TankUpdate tests ---

  #[test]
  fn test_update_tank_all_fields() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    let props_before = network.properties_version;

    network.update_tank("T1", &TankUpdate {
      elevation: Some(120.0),
      coordinates: Some((3.0, 4.0)),
      initial_level: Some(12.0),
      min_level: Some(2.0),
      max_level: Some(25.0),
      diameter: Some(60.0),
      min_volume: Some(5.0),
      overflow: Some(true),
    }).unwrap();

    let node = &network.nodes[0];
    assert_eq!(node.elevation, 120.0);
    assert_eq!(node.coordinates, Some((3.0, 4.0)));
    let NodeType::Tank(tank) = &node.node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.elevation, 120.0);
    assert_eq!(tank.initial_level, 12.0);
    assert_eq!(tank.min_level, 2.0);
    assert_eq!(tank.max_level, 25.0);
    assert_eq!(tank.diameter, 60.0);
    assert_eq!(tank.min_volume, 5.0);
    assert!(tank.overflow);
    assert_eq!(network.properties_version, props_before + 1);
    assert!(network.updated_nodes.contains(&0));
  }

  #[test]
  fn test_update_tank_partial_leaves_other_fields_unchanged() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();

    network.update_tank("T1", &TankUpdate {
      elevation: None,
      coordinates: None,
      initial_level: Some(15.0),
      min_level: None,
      max_level: None,
      diameter: None,
      min_volume: None,
      overflow: None,
    }).unwrap();

    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.initial_level, 15.0);
    assert_eq!(tank.min_level, 0.0);
    assert_eq!(tank.max_level, 20.0);
    assert_eq!(tank.diameter, 50.0);
    assert_eq!(network.nodes[0].elevation, 100.0);
  }

  #[test]
  fn test_update_tank_si_units() {
    // In CMH / SI, tank levels are stored internally in feet; user values should round-trip.
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    // add a tank directly in standard (ft) units so we have a known baseline
    network.add_node(test_tank_node("T1", 100.0 / MperFT)).unwrap();

    network.update_tank("T1", &TankUpdate {
      elevation: Some(50.0),
      coordinates: None,
      initial_level: Some(10.0),
      min_level: Some(0.0),
      max_level: Some(20.0),
      diameter: Some(20.0),
      min_volume: None,
      overflow: None,
    }).unwrap();

    // user-provided meters should be stored as feet internally
    assert!((network.nodes[0].elevation - 50.0 / MperFT).abs() < 1e-9);
    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert!((tank.elevation - 50.0 / MperFT).abs() < 1e-9);
    assert!((tank.initial_level - 10.0 / MperFT).abs() < 1e-9);
    assert!((tank.diameter - 20.0 / MperFT).abs() < 1e-9);
  }

  #[test]
  fn test_update_tank_wrong_node_type() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 0.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();

    let err = network.update_tank("J1", &TankUpdate {
      elevation: None,
      coordinates: None,
      initial_level: Some(1.0),
      min_level: None,
      max_level: None,
      diameter: None,
      min_volume: None,
      overflow: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotATank { .. }));
  }

  #[test]
  fn test_update_tank_not_found() {
    let mut network = Network::default();
    let err = network.update_tank("missing", &TankUpdate {
      elevation: None,
      coordinates: None,
      initial_level: None,
      min_level: None,
      max_level: None,
      diameter: None,
      min_volume: None,
      overflow: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
  }

  // --- ReservoirUpdate tests ---

  #[test]
  fn test_update_reservoir_elevation_and_coordinates() {
    let mut network = Network::default();
    network.add_node(test_reservoir_node("R1", 50.0)).unwrap();
    let props_before = network.properties_version;

    network.update_reservoir("R1", &ReservoirUpdate {
      elevation: Some(75.0),
      coordinates: Some((10.0, 20.0)),
      head_pattern: None,
    }).unwrap();

    assert_eq!(network.nodes[0].elevation, 75.0);
    assert_eq!(network.nodes[0].coordinates, Some((10.0, 20.0)));
    assert_eq!(network.properties_version, props_before + 1);
    assert!(network.updated_nodes.contains(&0));
  }

  #[test]
  fn test_update_reservoir_sets_head_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0, 2.0] });
    network.pattern_map.insert("P1".into(), 0);
    network.add_node(test_reservoir_node("R1", 50.0)).unwrap();

    network.update_reservoir("R1", &ReservoirUpdate {
      elevation: None,
      coordinates: None,
      head_pattern: Some(Some("P1".into())),
    }).unwrap();

    let NodeType::Reservoir(reservoir) = &network.nodes[0].node_type else {
      panic!("Expected Reservoir node type");
    };
    assert_eq!(reservoir.head_pattern.as_deref(), Some("P1"));
    assert_eq!(reservoir.head_pattern_index, Some(0));
  }

  #[test]
  fn test_update_reservoir_clears_head_pattern() {
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0] });
    network.pattern_map.insert("P1".into(), 0);
    let mut node = test_reservoir_node("R1", 50.0);
    if let NodeType::Reservoir(reservoir) = &mut node.node_type {
      reservoir.head_pattern = Some("P1".into());
      reservoir.head_pattern_index = Some(0);
    }
    network.add_node(node).unwrap();

    network.update_reservoir("R1", &ReservoirUpdate {
      elevation: None,
      coordinates: None,
      head_pattern: Some(None),
    }).unwrap();

    let NodeType::Reservoir(reservoir) = &network.nodes[0].node_type else {
      panic!("Expected Reservoir node type");
    };
    assert!(reservoir.head_pattern.is_none());
    assert!(reservoir.head_pattern_index.is_none());
  }

  #[test]
  fn test_update_reservoir_pattern_not_found() {
    let mut network = Network::default();
    network.add_node(test_reservoir_node("R1", 50.0)).unwrap();

    let err = network.update_reservoir("R1", &ReservoirUpdate {
      elevation: None,
      coordinates: None,
      head_pattern: Some(Some("missing".into())),
    }).unwrap_err();
    assert!(matches!(err, InputError::PatternNotFound { .. }));
  }

  #[test]
  fn test_update_reservoir_none_head_pattern_leaves_existing_pattern() {
    // `head_pattern: None` must be a no-op, not a clear.
    let mut network = Network::default();
    network.patterns.push(Pattern { id: "P1".into(), multipliers: vec![1.0] });
    network.pattern_map.insert("P1".into(), 0);
    let mut node = test_reservoir_node("R1", 50.0);
    if let NodeType::Reservoir(reservoir) = &mut node.node_type {
      reservoir.head_pattern = Some("P1".into());
      reservoir.head_pattern_index = Some(0);
    }
    network.add_node(node).unwrap();

    network.update_reservoir("R1", &ReservoirUpdate {
      elevation: Some(75.0),
      coordinates: None,
      head_pattern: None,
    }).unwrap();

    let NodeType::Reservoir(reservoir) = &network.nodes[0].node_type else {
      panic!("Expected Reservoir node type");
    };
    assert_eq!(reservoir.head_pattern.as_deref(), Some("P1"));
    assert_eq!(reservoir.head_pattern_index, Some(0));
    assert!((network.nodes[0].elevation - 75.0).abs() < 1e-9);
  }

  #[test]
  fn test_update_reservoir_wrong_node_type() {
    let mut network = Network::default();
    let data = JunctionData {
      elevation: 0.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    };
    network.add_junction("J1", &data).unwrap();

    let err = network.update_reservoir("J1", &ReservoirUpdate {
      elevation: None,
      coordinates: None,
      head_pattern: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotAReservoir { .. }));
  }

  #[test]
  fn test_update_reservoir_not_found() {
    let mut network = Network::default();
    let err = network.update_reservoir("missing", &ReservoirUpdate {
      elevation: None,
      coordinates: None,
      head_pattern: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
  }

  #[test]
  fn test_reset_changes_clears_updated_sets() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    network.update_tank("T1", &TankUpdate {
      elevation: None,
      coordinates: None,
      initial_level: Some(5.0),
      min_level: None,
      max_level: None,
      diameter: None,
      min_volume: None,
      overflow: None,
    }).unwrap();
    assert!(!network.updated_nodes.is_empty());

    network.reset_changes();
    assert!(network.updated_nodes.is_empty());
    assert!(network.updated_links.is_empty());
  }

  // --- Link helpers ---

  fn add_two_junctions(network: &mut Network) {
    for id in ["J1", "J2"] {
      network.add_junction(id, &JunctionData {
        elevation: 100.0,
        basedemand: 0.0,
        emitter_coefficient: 0.0,
        pattern: None,
        coordinates: None,
      }).unwrap();
    }
  }

  fn sample_pipe_data() -> PipeData {
    PipeData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      length: 1000.0,
      diameter: 12.0,
      roughness: 100.0,
      minor_loss: 0.0,
      check_valve: false,
      initial_status: LinkStatus::Open,
      vertices: None,
    }
  }

  // --- add_pipe tests ---

  #[test]
  fn test_add_pipe_basic() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    let topo_before = network.topology_version;

    network.add_pipe("P1", &sample_pipe_data()).unwrap();

    assert_eq!(network.links.len(), 1);
    assert_eq!(network.links[0].id, "P1".into());
    assert_eq!(network.links[0].start_node, 0);
    assert_eq!(network.links[0].end_node, 1);
    assert_eq!(network.link_map.get("P1").copied(), Some(0));
    // adding a link should bump the topology version
    assert_eq!(network.topology_version, topo_before + 1);
    let LinkType::Pipe(pipe) = &network.links[0].link_type else {
      panic!("Expected Pipe link type");
    };
    // diameter supplied in inches, stored internally in feet
    assert!((pipe.diameter - 1.0).abs() < 1e-9);
    assert_eq!(pipe.length, 1000.0);
    assert_eq!(pipe.headloss_formula, HeadlossFormula::HazenWilliams);
  }

  #[test]
  fn test_add_pipe_normalizes_minor_loss() {
    let mut network = Network::default();
    add_two_junctions(&mut network);

    network.add_pipe("P1", &PipeData {
      minor_loss: 2.0,
      ..sample_pipe_data()
    }).unwrap();

    let LinkType::Pipe(pipe) = &network.links[0].link_type else {
      panic!("Expected Pipe link type");
    };
    // internal representation: 0.02517 * K / d_ft^4 (d = 12in = 1ft)
    assert!((pipe.minor_loss - 0.05034).abs() < 1e-6);
  }

  #[test]
  fn test_add_pipe_si_units() {
    let mut network = Network::new(FlowUnits::CMH, HeadlossFormula::HazenWilliams);
    add_two_junctions(&mut network);

    network.add_pipe("P1", &PipeData {
      length: 1000.0,
      diameter: 300.0, // mm
      ..sample_pipe_data()
    }).unwrap();

    let LinkType::Pipe(pipe) = &network.links[0].link_type else {
      panic!("Expected Pipe link type");
    };
    // 300 mm -> 0.3 m -> ft
    assert!((pipe.diameter - (0.3 / MperFT)).abs() < 1e-9);
    assert!((pipe.length - (1000.0 / MperFT)).abs() < 1e-9);
  }

  #[test]
  fn test_add_pipe_unknown_node() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    let err = network.add_pipe("P1", &PipeData {
      start_node: "missing".into(),
      ..sample_pipe_data()
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
    assert_eq!(network.links.len(), 0);
  }

  #[test]
  fn test_add_pipe_duplicate_id() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let err = network.add_pipe("P1", &sample_pipe_data()).unwrap_err();
    assert!(matches!(err, InputError::LinkExists { .. }));
  }

  #[test]
  fn test_add_pipe_updates_tank_link_lists() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    add_two_junctions(&mut network);

    network.add_pipe("P1", &PipeData {
      start_node: "J1".into(),
      end_node: "T1".into(),
      ..sample_pipe_data()
    }).unwrap();

    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert_eq!(tank.links_to, vec![0]);
    assert!(tank.links_from.is_empty());
  }

  // --- add_pump tests ---

  fn push_head_curve(network: &mut Network) {
    network.curves.push(crate::model::curve::Curve {
      id: "C1".into(),
      x: vec![0.0, 500.0, 1000.0],
      y: vec![100.0, 80.0, 0.0],
    });
    network.curve_map.insert("C1".into(), 0);
  }

  #[test]
  fn test_add_pump_constant_power() {
    let mut network = Network::default();
    add_two_junctions(&mut network);

    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: None,
      power: 100.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.speed, 1.0);
    assert_eq!(pump.power, 100.0);
    assert!(pump.head_curve.is_none());
  }

  #[test]
  fn test_add_pump_with_head_curve() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);

    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("C1".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.head_curve_id.as_deref(), Some("C1"));
    assert!(pump.head_curve.is_some());
  }

  #[test]
  fn test_add_pump_unknown_curve() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    let err = network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("missing".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::CurveNotFound { .. }));
    assert_eq!(network.links.len(), 0);
  }

  // --- add_valve tests ---

  #[test]
  fn test_add_prv_applies_end_elevation_offset() {
    // user adds a PRV with setting 50 psi (US units). Internal setting should
    // be the converted pressure plus the end node's elevation (in feet).
    let mut network = Network::default();
    // J1 at 100 ft, J2 at 200 ft
    network.add_junction("J1", &JunctionData {
      elevation: 100.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0,
      basedemand: 0.0,
      emitter_coefficient: 0.0,
      pattern: None,
      coordinates: None,
    }).unwrap();

    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    // 50 psi -> 50 / PSIperFT ft, plus end-node elevation (J2 = 200 ft)
    let expected = 50.0 / crate::constants::PSIperFT + 200.0;
    assert!((valve.setting - expected).abs() < 1e-6);
    assert!(network.contains_pressure_control_valve);
  }

  #[test]
  fn test_add_tcv_no_elevation_offset() {
    let mut network = Network::default();
    add_two_junctions(&mut network);

    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::TCV,
      setting: 5.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.setting, 5.0);
    assert!(!network.contains_pressure_control_valve);
  }

  // --- update_pipe tests ---

  #[test]
  fn test_update_pipe_fields() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let props_before = network.properties_version;

    network.update_pipe("P1", &PipeUpdate {
      length: Some(500.0),
      diameter: Some(24.0), // inches -> 2 ft
      roughness: Some(120.0),
      minor_loss: Some(3.0),
      check_valve: Some(true),
    }).unwrap();

    let LinkType::Pipe(pipe) = &network.links[0].link_type else {
      panic!("Expected Pipe link type");
    };
    assert_eq!(pipe.length, 500.0);
    assert!((pipe.diameter - 2.0).abs() < 1e-9);
    assert_eq!(pipe.roughness, 120.0);
    assert!(pipe.check_valve);
    // minor_loss = 0.02517 * 3 / 2^4
    assert!((pipe.minor_loss - (0.02517 * 3.0 / 16.0)).abs() < 1e-9);
    assert_eq!(network.properties_version, props_before + 1);
    assert!(network.updated_links.contains(&0));
  }

  #[test]
  fn test_update_pipe_partial_preserves_minor_loss_coefficient() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &PipeData {
      minor_loss: 2.0,
      ..sample_pipe_data()
    }).unwrap();
    let k_expected = 0.02517 * 2.0 / 1.0_f64.powi(4);

    // change only the length; K should remain equivalent to 2.0 despite the
    // round-trip through convert_from_standard/convert_to_standard
    network.update_pipe("P1", &PipeUpdate {
      length: Some(250.0),
      diameter: None,
      roughness: None,
      minor_loss: None,
      check_valve: None,
    }).unwrap();

    let LinkType::Pipe(pipe) = &network.links[0].link_type else {
      panic!("Expected Pipe link type");
    };
    assert!((pipe.minor_loss - k_expected).abs() < 1e-9);
    assert_eq!(pipe.length, 250.0);
  }

  #[test]
  fn test_update_pipe_wrong_type() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: None,
      power: 10.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    let err = network.update_pipe("PU1", &PipeUpdate {
      length: Some(10.0),
      diameter: None,
      roughness: None,
      minor_loss: None,
      check_valve: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::LinkNotAPipe { .. }));
  }

  #[test]
  fn test_update_pipe_not_found() {
    let mut network = Network::default();
    let err = network.update_pipe("missing", &PipeUpdate {
      length: None, diameter: None, roughness: None,
      minor_loss: None, check_valve: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::LinkNotFound { .. }));
  }

  // --- update_pump tests ---

  #[test]
  fn test_update_pump_speed_and_power() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: None,
      power: 100.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    network.update_pump("PU1", &PumpUpdate {
      speed: Some(1.5),
      power: Some(200.0),
      head_curve_id: None,
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.speed, 1.5);
    assert_eq!(pump.power, 200.0);
    assert!(pump.head_curve_id.is_none());
    assert!(pump.head_curve.is_none());
  }

  #[test]
  fn test_update_pump_sets_head_curve() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: None,
      power: 10.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    network.update_pump("PU1", &PumpUpdate {
      speed: None,
      power: None,
      head_curve_id: Some(Some("C1".into())),
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.head_curve_id.as_deref(), Some("C1"));
    assert!(pump.head_curve.is_some());
  }

  #[test]
  fn test_update_pump_wrong_type() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();

    let err = network.update_pump("P1", &PumpUpdate {
      speed: Some(1.0), power: None, head_curve_id: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::LinkNotAPump { .. }));
  }

  #[test]
  fn test_update_pump_none_head_curve_leaves_existing_curve() {
    // `head_curve_id: None` must be a no-op, not a clear.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("C1".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    network.update_pump("PU1", &PumpUpdate {
      speed: Some(2.0),
      power: None,
      head_curve_id: None,
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert_eq!(pump.head_curve_id.as_deref(), Some("C1"));
    assert!(pump.head_curve.is_some());
    assert!((pump.speed - 2.0).abs() < 1e-9);
  }

  #[test]
  fn test_update_pump_clears_head_curve_explicitly() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("C1".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    network.update_pump("PU1", &PumpUpdate {
      speed: None,
      power: Some(50.0),
      head_curve_id: Some(None),
    }).unwrap();

    let LinkType::Pump(pump) = &network.links[0].link_type else {
      panic!("Expected Pump link type");
    };
    assert!(pump.head_curve_id.is_none());
    assert!(pump.head_curve.is_none());
  }

  // --- update_valve tests ---

  #[test]
  fn test_update_valve_setting_reapplies_prv_offset() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: Some(80.0),
      minor_loss: None,
      curve_id: None,
      valve_type: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    let expected = 80.0 / crate::constants::PSIperFT + 200.0;
    assert!((valve.setting - expected).abs() < 1e-6);
  }

  #[test]
  fn test_update_valve_diameter_preserves_prv_setting_semantics() {
    // If only the diameter is updated, the PRV setting (user-facing pressure)
    // should round-trip unchanged.
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();
    let setting_before = if let LinkType::Valve(valve) = &network.links[0].link_type {
      valve.setting
    } else { unreachable!() };

    network.update_valve("V1", &ValveUpdate {
      diameter: Some(24.0),
      setting: None,
      minor_loss: None,
      curve_id: None,
      valve_type: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert!((valve.setting - setting_before).abs() < 1e-9);
    assert!((valve.diameter - 2.0).abs() < 1e-9);
  }

  #[test]
  fn test_update_valve_wrong_type() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let err = network.update_valve("P1", &ValveUpdate {
      diameter: None, setting: None, minor_loss: None, curve_id: None,
      valve_type: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::LinkNotAValve { .. }));
  }

  #[test]
  fn test_update_valve_type_change_requires_setting() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::TCV,
      setting: 5.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    let err = network.update_valve("V1", &ValveUpdate {
      diameter: None, setting: None, minor_loss: None, curve_id: None,
      valve_type: Some(ValveType::PRV),
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_update_valve_type_change_to_gpv_requires_curve() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::TCV,
      setting: 5.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    let err = network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: Some(0.0),
      minor_loss: None,
      curve_id: None,
      valve_type: Some(ValveType::GPV),
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_update_valve_gpv_cannot_clear_curve() {
    // `Some(None)` on a GPV valve must be rejected — GPVs always need a curve.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::GPV,
      setting: 0.0,
      curve_id: Some("C1".into()),
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    let err = network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: None,
      minor_loss: None,
      curve_id: Some(None),
      valve_type: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_update_valve_gpv_none_curve_leaves_existing_curve() {
    // `curve_id: None` on a GPV valve without changing the type keeps the curve.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    push_head_curve(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::GPV,
      setting: 0.0,
      curve_id: Some("C1".into()),
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    network.update_valve("V1", &ValveUpdate {
      diameter: Some(10.0),
      setting: None,
      minor_loss: None,
      curve_id: None,
      valve_type: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.curve_id.as_deref(), Some("C1"));
    assert!(valve.valve_curve.is_some());
  }

  #[test]
  fn test_update_valve_tcv_none_curve_leaves_curve_none() {
    // Non-GPV case: `curve_id: None` is a no-op even when no curve is set.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::TCV,
      setting: 5.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();

    network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: Some(10.0),
      minor_loss: None,
      curve_id: None,
      valve_type: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert!(valve.curve_id.is_none());
    assert!(valve.valve_curve.is_none());
  }

  #[test]
  fn test_update_valve_type_change_tcv_to_prv_applies_offset() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::TCV,
      setting: 5.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();
    assert!(!network.contains_pressure_control_valve);

    network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: Some(50.0),
      minor_loss: None,
      curve_id: None,
      valve_type: Some(ValveType::PRV),
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.valve_type, ValveType::PRV);
    // 50 psi -> ft, plus end-node elevation (J2 = 200 ft)
    let expected = 50.0 / crate::constants::PSIperFT + 200.0;
    assert!((valve.setting - expected).abs() < 1e-6);
    assert!(network.contains_pressure_control_valve);
  }

  #[test]
  fn test_update_valve_type_change_prv_to_tcv_clears_flag_and_offset() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();
    assert!(network.contains_pressure_control_valve);

    network.update_valve("V1", &ValveUpdate {
      diameter: None,
      setting: Some(5.0),
      minor_loss: None,
      curve_id: None,
      valve_type: Some(ValveType::TCV),
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    assert_eq!(valve.valve_type, ValveType::TCV);
    // TCV settings are stored verbatim (no PSI->ft or elevation offset)
    assert_eq!(valve.setting, 5.0);
    assert!(!network.contains_pressure_control_valve);
  }

  // --- update_link tests ---

  #[test]
  fn test_update_link_vertices_and_status() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let topo_before = network.topology_version;
    let props_before = network.properties_version;

    network.update_link("P1", &LinkUpdate {
      start_node: None,
      end_node: None,
      vertices: Some(vec![(1.0, 2.0), (3.0, 4.0)]),
      initial_status: Some(LinkStatus::Closed),
    }).unwrap();

    assert_eq!(network.links[0].vertices.as_deref().unwrap(), &[(1.0, 2.0), (3.0, 4.0)]);
    assert_eq!(network.links[0].initial_status, LinkStatus::Closed);
    // vertices alone are visual; status bumps properties_version
    assert_eq!(network.topology_version, topo_before);
    assert_eq!(network.properties_version, props_before + 1);
    assert!(network.updated_links.contains(&0));
  }

  #[test]
  fn test_update_link_swap_endpoints_bumps_topology() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_junction("J3", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let topo_before = network.topology_version;

    network.update_link("P1", &LinkUpdate {
      start_node: None,
      end_node: Some("J3".into()),
      vertices: None,
      initial_status: None,
    }).unwrap();

    assert_eq!(network.links[0].end_node, 2);
    assert_eq!(network.links[0].end_node_id, "J3".into());
    assert_eq!(network.topology_version, topo_before + 1);
  }

  #[test]
  fn test_update_link_moves_tank_link_list_entry() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    network.add_node(test_tank_node("T2", 100.0)).unwrap();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();

    network.add_pipe("P1", &PipeData {
      start_node: "J1".into(),
      end_node: "T1".into(),
      ..sample_pipe_data()
    }).unwrap();

    network.update_link("P1", &LinkUpdate {
      start_node: None,
      end_node: Some("T2".into()),
      vertices: None,
      initial_status: None,
    }).unwrap();

    let NodeType::Tank(t1) = &network.nodes[0].node_type else { panic!(); };
    let NodeType::Tank(t2) = &network.nodes[1].node_type else { panic!(); };
    assert!(t1.links_to.is_empty(), "T1 should no longer be a destination");
    assert_eq!(t2.links_to, vec![0], "T2 should now be a destination");
  }

  #[test]
  fn test_update_link_prv_setting_follows_end_node_elevation() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 100.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 200.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J3", &JunctionData {
      elevation: 300.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();
    let setting_before = if let LinkType::Valve(valve) = &network.links[0].link_type {
      valve.setting
    } else { unreachable!() };

    network.update_link("V1", &LinkUpdate {
      start_node: None,
      end_node: Some("J3".into()),
      vertices: None,
      initial_status: None,
    }).unwrap();

    let LinkType::Valve(valve) = &network.links[0].link_type else {
      panic!("Expected Valve link type");
    };
    // the internal setting should have shifted by +100 ft (J3 - J2 elevations)
    assert!((valve.setting - (setting_before + 100.0)).abs() < 1e-6);
  }

  #[test]
  fn test_update_link_unknown_node() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let err = network.update_link("P1", &LinkUpdate {
      start_node: Some("missing".into()),
      end_node: None,
      vertices: None,
      initial_status: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
  }

  #[test]
  fn test_update_link_not_found() {
    let mut network = Network::default();
    let err = network.update_link("missing", &LinkUpdate {
      start_node: None, end_node: None, vertices: None, initial_status: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::LinkNotFound { .. }));
  }

  // --- remove_link tests ---

  #[test]
  fn test_remove_link_basic() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    let topo_before = network.topology_version;

    network.remove_link("P1", false).unwrap();

    assert_eq!(network.links.len(), 0);
    assert!(network.link_map.get("P1").is_none());
    assert_eq!(network.topology_version, topo_before + 1);
  }

  #[test]
  fn test_remove_link_reindexes_link_map_after_swap() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    network.add_pipe("P2", &sample_pipe_data()).unwrap();
    network.add_pipe("P3", &sample_pipe_data()).unwrap();

    // remove the middle link; swap_remove puts P3 into slot 1
    network.remove_link("P2", false).unwrap();

    assert_eq!(network.links.len(), 2);
    assert_eq!(network.link_map.get("P1").copied(), Some(0));
    assert_eq!(network.link_map.get("P3").copied(), Some(1));
    assert!(network.link_map.get("P2").is_none());
    assert_eq!(network.links[1].id, "P3".into());
  }

  #[test]
  fn test_remove_link_detaches_tank_links() {
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    add_two_junctions(&mut network);

    network.add_pipe("P1", &PipeData {
      start_node: "J1".into(),
      end_node: "T1".into(),
      ..sample_pipe_data()
    }).unwrap();

    network.remove_link("P1", false).unwrap();

    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    assert!(tank.links_to.is_empty());
    assert!(tank.links_from.is_empty());
  }

  #[test]
  fn test_remove_link_retargets_moved_link_in_tank_lists() {
    // Two tank-terminating pipes; remove the first so swap_remove moves
    // the last one into slot 0 and the tank lists must be retargeted.
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    add_two_junctions(&mut network);

    network.add_pipe("P1", &PipeData {
      start_node: "J1".into(),
      end_node: "T1".into(),
      ..sample_pipe_data()
    }).unwrap();
    network.add_pipe("P2", &PipeData {
      start_node: "J2".into(),
      end_node: "T1".into(),
      ..sample_pipe_data()
    }).unwrap();

    network.remove_link("P1", false).unwrap();

    let NodeType::Tank(tank) = &network.nodes[0].node_type else {
      panic!("Expected Tank node type");
    };
    // P2 is now at index 0
    assert_eq!(tank.links_to, vec![0]);
    assert_eq!(network.links[0].id, "P2".into());
  }

  #[test]
  fn test_remove_link_clears_pressure_control_flag() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_valve("V1", &ValveData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      diameter: 12.0,
      valve_type: ValveType::PRV,
      setting: 50.0,
      curve_id: None,
      minor_loss: 0.0,
      initial_status: LinkStatus::Active,
      vertices: None,
    }).unwrap();
    assert!(network.contains_pressure_control_valve);

    network.remove_link("V1", false).unwrap();
    assert!(!network.contains_pressure_control_valve);
  }

  #[test]
  fn test_remove_link_referenced_by_control_fails() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    network.controls.push(Control {
      condition: ControlCondition::Time { seconds: 0 },
      link_id: "P1".into(),
      setting: None,
      status: Some(LinkStatus::Closed),
    });

    let err = network.remove_link("P1", false).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
    assert_eq!(network.links.len(), 1);
  }

  #[test]
  fn test_remove_link_not_found() {
    let mut network = Network::default();
    let err = network.remove_link("missing", false).unwrap_err();
    assert!(matches!(err, InputError::LinkNotFound { .. }));
  }

  // --- remove_node tests ---

  #[test]
  fn test_remove_node_basic() {
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    let topo_before = network.topology_version;

    network.remove_node("J1", false).unwrap();

    assert_eq!(network.nodes.len(), 0);
    assert!(network.node_map.get("J1").is_none());
    assert_eq!(network.topology_version, topo_before + 1);
  }

  #[test]
  fn test_remove_node_reindexes_links_after_swap() {
    // J1, J2, J3 with a pipe J1->J3. Remove J2; J3 is swapped into slot 1
    // and the pipe's end_node index must follow.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_junction("J3", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_pipe("P1", &PipeData {
      start_node: "J1".into(),
      end_node: "J3".into(),
      ..sample_pipe_data()
    }).unwrap();

    network.remove_node("J2", false).unwrap();

    assert_eq!(network.nodes.len(), 2);
    assert_eq!(network.node_map.get("J1").copied(), Some(0));
    assert_eq!(network.node_map.get("J3").copied(), Some(1));
    assert_eq!(network.links[0].start_node, 0);
    assert_eq!(network.links[0].end_node, 1);
  }

  #[test]
  fn test_remove_node_cascade_removes_connected_links() {
    // J1 -P1- J2 -P2- J3 and a stray P3 (J1 -> J3). Removing J1 should
    // cascade-remove P1 and P3 but leave P2 alone.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_junction("J3", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    network.add_pipe("P2", &PipeData {
      start_node: "J2".into(), end_node: "J3".into(),
      ..sample_pipe_data()
    }).unwrap();
    network.add_pipe("P3", &PipeData {
      start_node: "J1".into(), end_node: "J3".into(),
      ..sample_pipe_data()
    }).unwrap();

    network.remove_node("J1", false).unwrap();

    assert_eq!(network.nodes.len(), 2);
    assert!(network.node_map.get("J1").is_none());
    assert_eq!(network.links.len(), 1);
    assert!(network.link_map.get("P1").is_none());
    assert!(network.link_map.get("P3").is_none());
    assert!(network.link_map.get("P2").is_some());
  }

  #[test]
  fn test_remove_node_cascade_fails_if_link_referenced_by_control() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_pipe("P1", &sample_pipe_data()).unwrap();
    network.controls.push(Control {
      condition: ControlCondition::Time { seconds: 0 },
      link_id: "P1".into(),
      setting: None,
      status: Some(LinkStatus::Closed),
    });
    let err = network.remove_node("J1", false).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_remove_node_referenced_by_control_fails() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.controls.push(Control {
      condition: ControlCondition::HighPressure { node_index: 0, target: 50.0 },
      link_id: "dummy".into(),
      setting: None,
      status: Some(LinkStatus::Closed),
    });

    let err = network.remove_node("J1", false).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_remove_node_reindexes_control_node_index() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    network.add_junction("J2", &JunctionData {
      elevation: 0.0, basedemand: 0.0, emitter_coefficient: 0.0,
      pattern: None, coordinates: None,
    }).unwrap();
    // control targets J2 (index 1)
    network.controls.push(Control {
      condition: ControlCondition::HighPressure { node_index: 1, target: 50.0 },
      link_id: "dummy".into(),
      setting: None,
      status: Some(LinkStatus::Closed),
    });

    // remove J1 (index 0); J2 gets swapped into slot 0, control must follow
    network.remove_node("J1", false).unwrap();

    let ControlCondition::HighPressure { node_index, .. } = network.controls[0].condition else {
      panic!("expected HighPressure control");
    };
    assert_eq!(node_index, 0);
  }

  #[test]
  fn test_remove_node_unconditional_drops_referencing_controls() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    add_two_junctions(&mut network);
    // unrelated control on J2 (index 1) must survive and get reindexed to 0
    network.controls.push(Control {
      condition: ControlCondition::HighPressure { node_index: 1, target: 50.0 },
      link_id: "dummy".into(),
      setting: None,
      status: Some(LinkStatus::Closed),
    });
    // two controls (HighPressure + LowPressure) targeting J1 must be dropped
    network.controls.push(Control {
      condition: ControlCondition::HighPressure { node_index: 0, target: 60.0 },
      link_id: "dummy".into(), setting: None, status: Some(LinkStatus::Closed),
    });
    network.controls.push(Control {
      condition: ControlCondition::LowPressure { node_index: 0, target: 10.0 },
      link_id: "dummy".into(), setting: None, status: Some(LinkStatus::Open),
    });

    network.remove_node("J1", true).unwrap();

    // only the J2-targeting control remains; its node_index is reindexed
    // from 1 to 0 because swap_remove moved J2 into J1's old slot.
    assert_eq!(network.controls.len(), 1);
    let ControlCondition::HighPressure { node_index, .. } = network.controls[0].condition else {
      panic!("expected HighPressure control");
    };
    assert_eq!(node_index, 0);
  }

  #[test]
  fn test_remove_node_unconditional_drops_tank_level_controls() {
    use crate::model::control::{Control, ControlCondition};
    let mut network = Network::default();
    network.add_node(test_tank_node("T1", 100.0)).unwrap();
    network.controls.push(Control {
      condition: ControlCondition::HighLevel { tank_index: 0, target: 20.0 },
      link_id: "dummy".into(), setting: None, status: Some(LinkStatus::Closed),
    });
    network.controls.push(Control {
      condition: ControlCondition::LowLevel { tank_index: 0, target: 2.0 },
      link_id: "dummy".into(), setting: None, status: Some(LinkStatus::Open),
    });

    network.remove_node("T1", true).unwrap();

    assert!(network.controls.is_empty());
    assert!(network.node_map.get("T1").is_none());
  }

  #[test]
  fn test_remove_node_not_found() {
    let mut network = Network::default();
    let err = network.remove_node("missing", false).unwrap_err();
    assert!(matches!(err, InputError::NodeNotFound { .. }));
  }

  // --- pattern tests ---

  #[test]
  fn test_add_pattern_basic() {
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData {
      multipliers: vec![1.0, 1.2, 0.8],
    }).unwrap();
    assert_eq!(network.patterns.len(), 1);
    assert_eq!(network.pattern_map.get("P1").copied(), Some(0));
    assert_eq!(network.patterns[0].multipliers, vec![1.0, 1.2, 0.8]);
  }

  #[test]
  fn test_add_pattern_duplicate_fails() {
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData { multipliers: vec![1.0] }).unwrap();
    let err = network.add_pattern("P1", &PatternData { multipliers: vec![2.0] }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_update_pattern_replaces_multipliers() {
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData { multipliers: vec![1.0] }).unwrap();
    let props_before = network.properties_version;

    network.update_pattern("P1", &PatternUpdate {
      multipliers: Some(vec![0.5, 1.5]),
    }).unwrap();

    assert_eq!(network.patterns[0].multipliers, vec![0.5, 1.5]);
    assert_eq!(network.properties_version, props_before + 1);
  }

  #[test]
  fn test_update_pattern_not_found() {
    let mut network = Network::default();
    let err = network.update_pattern("missing", &PatternUpdate::default()).unwrap_err();
    assert!(matches!(err, InputError::PatternNotFound { .. }));
  }

  #[test]
  fn test_remove_pattern_basic() {
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData { multipliers: vec![1.0] }).unwrap();
    network.remove_pattern("P1").unwrap();
    assert!(network.patterns.is_empty());
    assert!(network.pattern_map.get("P1").is_none());
  }

  #[test]
  fn test_remove_pattern_reindexes_junctions_and_reservoirs() {
    // two patterns P1, P2. Junction J1 uses P2; remove P1 and verify P2
    // now sits at index 0 AND J1's cached pattern_index points to 0.
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData { multipliers: vec![1.0] }).unwrap();
    network.add_pattern("P2", &PatternData { multipliers: vec![2.0] }).unwrap();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 1.0, emitter_coefficient: 0.0,
      pattern: Some("P2".into()),
      coordinates: None,
    }).unwrap();
    network.add_reservoir("R1", &ReservoirData {
      elevation: 100.0,
      head_pattern: Some("P2".into()),
      coordinates: None,
    }).unwrap();

    network.remove_pattern("P1").unwrap();

    assert_eq!(network.pattern_map.get("P2").copied(), Some(0));
    let NodeType::Junction(j) = &network.nodes[0].node_type else { panic!() };
    assert_eq!(j.pattern_index, Some(0));
    let NodeType::Reservoir(r) = &network.nodes[1].node_type else { panic!() };
    assert_eq!(r.head_pattern_index, Some(0));
  }

  #[test]
  fn test_remove_pattern_referenced_fails() {
    let mut network = Network::default();
    network.add_pattern("P1", &PatternData { multipliers: vec![1.0] }).unwrap();
    network.add_junction("J1", &JunctionData {
      elevation: 0.0, basedemand: 1.0, emitter_coefficient: 0.0,
      pattern: Some("P1".into()),
      coordinates: None,
    }).unwrap();
    let err = network.remove_pattern("P1").unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
    assert_eq!(network.patterns.len(), 1);
  }

  #[test]
  fn test_remove_pattern_not_found() {
    let mut network = Network::default();
    let err = network.remove_pattern("missing").unwrap_err();
    assert!(matches!(err, InputError::PatternNotFound { .. }));
  }

  // --- curve tests ---

  #[test]
  fn test_add_curve_basic() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 10.0, 20.0],
      y: vec![100.0, 80.0, 40.0],
    }).unwrap();
    assert_eq!(network.curves.len(), 1);
    assert_eq!(network.curve_map.get("C1").copied(), Some(0));
  }

  #[test]
  fn test_add_curve_rejects_mismatched_lengths() {
    let mut network = Network::default();
    let err = network.add_curve("C1", &CurveData {
      x: vec![0.0, 10.0],
      y: vec![100.0],
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_add_curve_rejects_non_monotonic_x() {
    let mut network = Network::default();
    let err = network.add_curve("C1", &CurveData {
      x: vec![0.0, 10.0, 5.0],
      y: vec![1.0, 2.0, 3.0],
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_add_curve_rejects_empty() {
    let mut network = Network::default();
    let err = network.add_curve("C1", &CurveData {
      x: vec![], y: vec![],
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_add_curve_duplicate_fails() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0], y: vec![1.0, 0.0],
    }).unwrap();
    let err = network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0], y: vec![1.0, 0.0],
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
  }

  #[test]
  fn test_update_curve_rebuilds_pump_head_curve() {
    // a pump referencing C1 should have its cached HeadCurve refreshed when
    // C1's axes are updated, and the owning link marked as updated.
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 50.0, 100.0],
      y: vec![100.0, 80.0, 40.0],
    }).unwrap();
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("C1".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();
    network.reset_changes();
    let props_before = network.properties_version;

    network.update_curve("C1", &CurveUpdate {
      x: Some(vec![0.0, 25.0, 60.0]),
      y: Some(vec![120.0, 90.0, 30.0]),
    }).unwrap();

    assert_eq!(network.curves[0].x, vec![0.0, 25.0, 60.0]);
    assert_eq!(network.curves[0].y, vec![120.0, 90.0, 30.0]);
    let LinkType::Pump(pump) = &network.links[0].link_type else { panic!() };
    assert!(pump.head_curve.is_some());
    // the refreshed shutoff head (in feet) should reflect the new y[0]
    let hc = pump.head_curve.as_ref().unwrap();
    assert!((hc.statistics.h_shutoff - 120.0).abs() < 1e-6);
    assert!(network.updated_links.contains(&0));
    assert_eq!(network.properties_version, props_before + 1);
  }

  #[test]
  fn test_update_curve_partial_preserves_other_axis() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0, 2.0],
      y: vec![10.0, 8.0, 5.0],
    }).unwrap();

    network.update_curve("C1", &CurveUpdate {
      x: None,
      y: Some(vec![20.0, 15.0, 9.0]),
    }).unwrap();

    assert_eq!(network.curves[0].x, vec![0.0, 1.0, 2.0]);
    assert_eq!(network.curves[0].y, vec![20.0, 15.0, 9.0]);
  }

  #[test]
  fn test_update_curve_rejects_invalid_combined_axes() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0, 2.0],
      y: vec![10.0, 8.0, 5.0],
    }).unwrap();

    let err = network.update_curve("C1", &CurveUpdate {
      x: Some(vec![0.0, 2.0]), // shorter than existing y
      y: None,
    }).unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
    // curve must be unchanged on failure
    assert_eq!(network.curves[0].x, vec![0.0, 1.0, 2.0]);
  }

  #[test]
  fn test_update_curve_not_found() {
    let mut network = Network::default();
    let err = network.update_curve("missing", &CurveUpdate::default()).unwrap_err();
    assert!(matches!(err, InputError::CurveNotFound { .. }));
  }

  #[test]
  fn test_remove_curve_basic() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0], y: vec![1.0, 0.0],
    }).unwrap();
    network.remove_curve("C1").unwrap();
    assert!(network.curves.is_empty());
    assert!(network.curve_map.get("C1").is_none());
  }

  #[test]
  fn test_remove_curve_reindexes_curve_map() {
    let mut network = Network::default();
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 1.0], y: vec![1.0, 0.0],
    }).unwrap();
    network.add_curve("C2", &CurveData {
      x: vec![0.0, 1.0], y: vec![2.0, 1.0],
    }).unwrap();

    network.remove_curve("C1").unwrap();
    assert_eq!(network.curves.len(), 1);
    assert_eq!(network.curve_map.get("C2").copied(), Some(0));
  }

  #[test]
  fn test_remove_curve_referenced_by_pump_fails() {
    let mut network = Network::default();
    add_two_junctions(&mut network);
    network.add_curve("C1", &CurveData {
      x: vec![0.0, 50.0, 100.0],
      y: vec![100.0, 80.0, 40.0],
    }).unwrap();
    network.add_pump("PU1", &PumpData {
      start_node: "J1".into(),
      end_node: "J2".into(),
      speed: 1.0,
      head_curve_id: Some("C1".into()),
      power: 0.0,
      initial_status: LinkStatus::Open,
      vertices: None,
    }).unwrap();

    let err = network.remove_curve("C1").unwrap_err();
    assert!(matches!(err, InputError::Parse { .. }));
    assert_eq!(network.curves.len(), 1);
  }

  #[test]
  fn test_remove_curve_not_found() {
    let mut network = Network::default();
    let err = network.remove_curve("missing").unwrap_err();
    assert!(matches!(err, InputError::CurveNotFound { .. }));
  }
}