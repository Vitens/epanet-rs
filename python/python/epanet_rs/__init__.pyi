EN_JUNCTION: int
EN_RESERVOIR: int
EN_TANK: int
EN_CVPIPE: int
EN_PIPE: int
EN_PUMP: int
EN_PRV: int
EN_PSV: int
EN_PBV: int
EN_FCV: int
EN_TCV: int
EN_GPV: int
EN_PCV: int
EN_NODECOUNT: int
EN_TANKCOUNT: int
EN_LINKCOUNT: int
EN_PATCOUNT: int
EN_CURVECOUNT: int
EN_CONTROLCOUNT: int
EN_RULECOUNT: int
EN_ELEVATION: int
EN_BASEDEMAND: int
EN_PATTERN: int
EN_EMITTER: int
EN_INITQUAL: int
EN_SOURCEQUAL: int
EN_SOURCEPAT: int
EN_SOURCETYPE: int
EN_TANKLEVEL: int
EN_DEMAND: int
EN_HEAD: int
EN_PRESSURE: int
EN_QUALITY: int
EN_SOURCEMASS: int
EN_INITVOLUME: int
EN_MIXMODEL: int
EN_MIXZONEVOL: int
EN_TANKDIAM: int
EN_MINVOLUME: int
EN_VOLCURVE: int
EN_MINLEVEL: int
EN_MAXLEVEL: int
EN_MIXFRACTION: int
EN_TANK_KBULK: int
EN_TANKVOLUME: int
EN_MAXVOLUME: int
EN_CANOVERFLOW: int
EN_DEMANDDEFICIT: int
EN_DIAMETER: int
EN_LENGTH: int
EN_ROUGHNESS: int
EN_MINORLOSS: int
EN_INITSTATUS: int
EN_INITSETTING: int
EN_KBULK: int
EN_KWALL: int
EN_FLOW: int
EN_VELOCITY: int
EN_HEADLOSS: int
EN_STATUS: int
EN_SETTING: int
EN_ENERGY: int
EN_LINKQUAL: int
EN_LINKPATTERN: int
EN_PUMPSTATE: int
EN_PUMPEFFIC: int
EN_PUMPPOWER: int
EN_PUMPHCURVE: int
EN_PUMPECURVE: int
EN_PUMPECOST: int
EN_PUMPEPAT: int
EN_DURATION: int
EN_HYDSTEP: int
EN_QUALSTEP: int
EN_PATTERNSTEP: int
EN_PATTERNSTART: int
EN_REPORTSTEP: int
EN_REPORTSTART: int
EN_RULESTEP: int
EN_STATISTIC: int
EN_PERIODS: int
EN_STARTTIME: int
EN_HTIME: int
EN_QTIME: int
EN_HALTFLAG: int
EN_NEXTEVENT: int
EN_NEXTEVENTTANK: int
EN_TRIALS: int
EN_ACCURACY: int
EN_TOLERANCE: int
EN_EMITEXPON: int
EN_DEMANDMULT: int
EN_HEADERROR: int
EN_FLOWCHANGE: int
EN_HEADLOSSFORM: int
EN_GLOBALEFFIC: int
EN_GLOBALPRICE: int
EN_GLOBALPATTERN: int
EN_DEMANDCHARGE: int
EN_SP_GRAVITY: int
EN_SP_VISCOS: int
EN_UNBALANCED: int
EN_CHECKFREQ: int
EN_MAXCHECK: int
EN_DAMPLIMIT: int
EN_SP_DIFFUS: int
EN_BULKORDER: int
EN_WALLORDER: int
EN_TANKORDER: int
EN_CONCENLIMIT: int
EN_DDA: int
EN_PDA: int
EN_HW: int
EN_DW: int
EN_CM: int
EN_CFS: int
EN_GPM: int
EN_MGD: int
EN_IMGD: int
EN_AFD: int
EN_LPS: int
EN_LPM: int
EN_MLD: int
EN_CMH: int
EN_CMD: int
EN_MISSING: float

class SolverResult:
    """Results from a hydraulic simulation."""

    @property
    def flows(self) -> list[list[float]]: ...
    @property
    def heads(self) -> list[list[float]]: ...
    @property
    def demands(self) -> list[list[float]]: ...

class Project:
    """EPANET project handle — mirrors the EPANET 2.3 toolkit API."""

    def __init__(self) -> None: ...
    def open(self, inp_file: str) -> None: ...
    def close(self) -> None: ...
    def saveinpfile(self, path: str) -> None: ...

    # Hydraulic solver
    def openH(self) -> None: ...
    def initH(self, initflag: int = 0) -> None: ...
    def runH(self, time: int) -> None: ...
    def nextH(self) -> int: ...
    def solveH(self, parallel: bool = False) -> SolverResult: ...
    def closeH(self) -> None: ...

    # Water quality (stubs)
    def openQ(self) -> None: ...
    def initQ(self, initflag: int = 0) -> None: ...
    def closeQ(self) -> None: ...

    # Counts
    def getcount(self, object: int) -> int: ...

    # Nodes
    def addnode(self, id: str, node_type: int) -> int: ...
    def getnodeindex(self, id: str) -> int: ...
    def getnodeid(self, index: int) -> str: ...
    def setnodeid(self, index: int, id: str) -> None: ...
    def getnodetype(self, index: int) -> int: ...
    def getnodevalue(self, index: int, property: int) -> float: ...
    def setnodevalue(self, index: int, property: int, value: float) -> None: ...
    def getcoord(self, index: int) -> tuple[float, float]: ...
    def setcoord(self, index: int, x: float, y: float) -> None: ...
    def deletenode(self, index: int, action_code: int = 0) -> None: ...

    # Links
    def addlink(self, id: str, link_type: int, start_node: str, end_node: str) -> int: ...
    def deletelink(self, index: int, action_code: int = 0) -> None: ...
    def getlinkindex(self, id: str) -> int: ...
    def getlinkid(self, index: int) -> str: ...
    def setlinkid(self, index: int, id: str) -> None: ...
    def getlinktype(self, index: int) -> int: ...
    def getlinknodes(self, index: int) -> tuple[int, int]: ...
    def setlinknodes(self, index: int, start_node: int, end_node: int) -> None: ...
    def getlinkvalue(self, index: int, property: int) -> float: ...
    def setlinkvalue(self, index: int, property: int, value: float) -> None: ...
    def getheadcurveindex(self, index: int) -> int: ...
    def setheadcurveindex(self, index: int, head_curve_index: int) -> None: ...

    # Patterns
    def addpattern(self, id: str) -> None: ...
    def deletepattern(self, index: int) -> None: ...
    def getpatternindex(self, id: str) -> int: ...
    def getpatternid(self, index: int) -> str: ...
    def setpatternid(self, index: int, id: str) -> None: ...
    def getpatternlen(self, index: int) -> int: ...
    def getpatternvalue(self, index: int, time: int) -> float: ...
    def setpattern(self, index: int, multipliers: list[float]) -> None: ...
    def getaveragepatternvalue(self, index: int) -> float: ...

    # Curves
    def addcurve(self, id: str) -> None: ...
    def getcurveindex(self, id: str) -> int: ...
    def getcurveid(self, index: int) -> str: ...
    def setcurveid(self, index: int, id: str) -> None: ...
    def getcurvelen(self, index: int) -> int: ...
    def getcurvevalue(self, index: int, point_index: int) -> tuple[float, float]: ...
    def getcurve(self, index: int) -> tuple[list[float], list[float]]: ...
    def setcurve(self, index: int, x: list[float], y: list[float]) -> None: ...

    # Time parameters
    def settimeparam(self, param: int, value: int) -> None: ...

    # Options
    def getoption(self, option: int) -> float: ...
    def setoption(self, option: int, value: float) -> None: ...
    def setdemandmodel(
        self,
        demand_model: int,
        minimum_pressure: float,
        required_pressure: float,
        pressure_exponent: float,
    ) -> None: ...

    # Error
    @staticmethod
    def geterror(errcode: int) -> str: ...
