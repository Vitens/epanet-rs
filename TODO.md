# Ideas and To-Do

- [x] Refactor Solver to output vector of heads/flows instead of in-place assignment (allow parallization of solve loop)
- [x] Parallelize sequential solves with Rayon if possible (no tanks)
- [x] Refactor network build from inp file, split logic into separate methods
- [x] Refactor network, link and node setup, split into separate files and compute coefficients using traits
- [x] Implement Darcy Weisbach
- [x] Input file support for PATTERN, TIMES, OPTIONS
- [x] Support for PATTERNS in solve loop 
- [x] Implement unit conversion
- [x] Implement pumps
- [x] Single point curves
- [x] Three point curves
- [x] Custom curves
- [x] Use serde for serializing to/from json and msgpack files
- [x] Implement FCV, TCV, FCV, PBV
- [x] Implement PRV and PSV valves
- [x] Coeff stability checks (pipe, valve, pump) (RQTol)
- [x] FCV with Curve
- [x] Implement GPV
- [x] Implement tank logic
- [x] Fixed power pumps
- [x] Implement CONTROLS
- [x] Refactor how time steps are calculated to match EPANET
- [x] Implement EMITTERS
- [x] Pressure Driven Analysis
- [x] Get rid of all the Panics and Unwraps, use Result instead
- [x] Network input file writing
- [x] Fix Unit conversion when serializing using Serde
- [x] Unit conversion for GPV head curves
- [x] Setup unit and integration testing, improve coverage
- [x] Setup automated benchmarks
- [ ] Create EPANET2_3 API compatibility layer
- [ ] Implement headerror convergence criteria
- [ ] Update convergence criteria for PDA
- [ ] Implement RULES
- [ ] Implement report start and statistic option
- [ ] Tanks with volume curve
- [ ] Network validation (connectivity, presence of tanks)
- [ ] Identify and remove unconnected zones from the solver

# API Methods
- [ ] fix pattern removal (should disasocciate nodes from pattern)

# Cleanup and refactoring
- [ ] Refactor curves to use a combined x and y vector instead of separate vectors
- [ ] Investigate performance degradation using HashMap for lookup of nodes, links, curves, patterns instead of Vec index
- [ ] Otherwise, write a method to remap all indices based on ids
- [ ] Setting tank_level should be possible?
- [ ] Implement get_* for getting link/node/curve/pattern/control/rule properties in correct units