# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.31
#
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _simV2
import new
new_instancemethod = new.instancemethod
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'PySwigObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


VisItDisconnect = _simV2.VisItDisconnect
VisIt_VariableData_setDataC = _simV2.VisIt_VariableData_setDataC
VisIt_VariableData_setDataI = _simV2.VisIt_VariableData_setDataI
VisIt_VariableData_setDataF = _simV2.VisIt_VariableData_setDataF
VisIt_VariableData_setDataD = _simV2.VisIt_VariableData_setDataD
VisItProcessEngineCommand = _simV2.VisItProcessEngineCommand
VisItSetBroadcastIntFunction = _simV2.VisItSetBroadcastIntFunction
VisItSetBroadcastStringFunction = _simV2.VisItSetBroadcastStringFunction
VisItSetParallel = _simV2.VisItSetParallel
VisItSetParallelRank = _simV2.VisItSetParallelRank
VisItSetDirectory = _simV2.VisItSetDirectory
VisItSetOptions = _simV2.VisItSetOptions
VisItSetupEnvironment = _simV2.VisItSetupEnvironment
VisItSetupEnvironment2 = _simV2.VisItSetupEnvironment2
VisItGetEnvironment = _simV2.VisItGetEnvironment
VisItInitializeSocketAndDumpSimFile = _simV2.VisItInitializeSocketAndDumpSimFile
VisItDetectInput = _simV2.VisItDetectInput
VisItDetectInputWithTimeout = _simV2.VisItDetectInputWithTimeout
VisItAttemptToCompleteConnection = _simV2.VisItAttemptToCompleteConnection
VisItSetSlaveProcessCallback = _simV2.VisItSetSlaveProcessCallback
VisItSetCommandCallback = _simV2.VisItSetCommandCallback
VisItTimeStepChanged = _simV2.VisItTimeStepChanged
VisItUpdatePlots = _simV2.VisItUpdatePlots
VisItExecuteCommand = _simV2.VisItExecuteCommand
VisItIsConnected = _simV2.VisItIsConnected
VisItGetLastError = _simV2.VisItGetLastError
VisItSynchronize = _simV2.VisItSynchronize
VisItEnableSynchronize = _simV2.VisItEnableSynchronize
VisItDebug1 = _simV2.VisItDebug1
VisItDebug2 = _simV2.VisItDebug2
VisItDebug3 = _simV2.VisItDebug3
VisItDebug4 = _simV2.VisItDebug4
VisItDebug5 = _simV2.VisItDebug5
VisItOpenTraceFile = _simV2.VisItOpenTraceFile
VisItCloseTraceFile = _simV2.VisItCloseTraceFile
VisItSaveWindow = _simV2.VisItSaveWindow
VisItSetActivateTimestep = _simV2.VisItSetActivateTimestep
VisItSetGetMetaData = _simV2.VisItSetGetMetaData
VisItSetGetMesh = _simV2.VisItSetGetMesh
VisItSetGetMaterial = _simV2.VisItSetGetMaterial
VisItSetGetSpecies = _simV2.VisItSetGetSpecies
VisItSetGetVariable = _simV2.VisItSetGetVariable
VisItSetGetMixedVariable = _simV2.VisItSetGetMixedVariable
VisItSetGetCurve = _simV2.VisItSetGetCurve
VisItSetGetDomainList = _simV2.VisItSetGetDomainList
VisItSetGetDomainBoundaries = _simV2.VisItSetGetDomainBoundaries
VisItSetGetDomainNesting = _simV2.VisItSetGetDomainNesting
VISIT_INVALID_HANDLE = _simV2.VISIT_INVALID_HANDLE
VISIT_ERROR = _simV2.VISIT_ERROR
VISIT_OKAY = _simV2.VISIT_OKAY
VISIT_NODATA = _simV2.VISIT_NODATA
VISIT_MESHTYPE_UNKNOWN = _simV2.VISIT_MESHTYPE_UNKNOWN
VISIT_MESHTYPE_RECTILINEAR = _simV2.VISIT_MESHTYPE_RECTILINEAR
VISIT_MESHTYPE_CURVILINEAR = _simV2.VISIT_MESHTYPE_CURVILINEAR
VISIT_MESHTYPE_UNSTRUCTURED = _simV2.VISIT_MESHTYPE_UNSTRUCTURED
VISIT_MESHTYPE_POINT = _simV2.VISIT_MESHTYPE_POINT
VISIT_MESHTYPE_CSG = _simV2.VISIT_MESHTYPE_CSG
VISIT_MESHTYPE_AMR = _simV2.VISIT_MESHTYPE_AMR
VISIT_VARCENTERING_NODE = _simV2.VISIT_VARCENTERING_NODE
VISIT_VARCENTERING_ZONE = _simV2.VISIT_VARCENTERING_ZONE
VISIT_VARTYPE_UNKNOWN = _simV2.VISIT_VARTYPE_UNKNOWN
VISIT_VARTYPE_SCALAR = _simV2.VISIT_VARTYPE_SCALAR
VISIT_VARTYPE_VECTOR = _simV2.VISIT_VARTYPE_VECTOR
VISIT_VARTYPE_TENSOR = _simV2.VISIT_VARTYPE_TENSOR
VISIT_VARTYPE_SYMMETRIC_TENSOR = _simV2.VISIT_VARTYPE_SYMMETRIC_TENSOR
VISIT_VARTYPE_MATERIAL = _simV2.VISIT_VARTYPE_MATERIAL
VISIT_VARTYPE_MATSPECIES = _simV2.VISIT_VARTYPE_MATSPECIES
VISIT_VARTYPE_LABEL = _simV2.VISIT_VARTYPE_LABEL
VISIT_VARTYPE_ARRAY = _simV2.VISIT_VARTYPE_ARRAY
VISIT_VARTYPE_MESH = _simV2.VISIT_VARTYPE_MESH
VISIT_VARTYPE_CURVE = _simV2.VISIT_VARTYPE_CURVE
VISIT_CMDARG_NONE = _simV2.VISIT_CMDARG_NONE
VISIT_CMDARG_INT = _simV2.VISIT_CMDARG_INT
VISIT_CMDARG_FLOAT = _simV2.VISIT_CMDARG_FLOAT
VISIT_CMDARG_STRING = _simV2.VISIT_CMDARG_STRING
VISIT_SIMMODE_UNKNOWN = _simV2.VISIT_SIMMODE_UNKNOWN
VISIT_SIMMODE_RUNNING = _simV2.VISIT_SIMMODE_RUNNING
VISIT_SIMMODE_STOPPED = _simV2.VISIT_SIMMODE_STOPPED
VISIT_DATATYPE_CHAR = _simV2.VISIT_DATATYPE_CHAR
VISIT_DATATYPE_INT = _simV2.VISIT_DATATYPE_INT
VISIT_DATATYPE_FLOAT = _simV2.VISIT_DATATYPE_FLOAT
VISIT_DATATYPE_DOUBLE = _simV2.VISIT_DATATYPE_DOUBLE
VISIT_OWNER_SIM = _simV2.VISIT_OWNER_SIM
VISIT_OWNER_VISIT = _simV2.VISIT_OWNER_VISIT
VISIT_OWNER_COPY = _simV2.VISIT_OWNER_COPY
VISIT_CELL_BEAM = _simV2.VISIT_CELL_BEAM
VISIT_CELL_TRI = _simV2.VISIT_CELL_TRI
VISIT_CELL_QUAD = _simV2.VISIT_CELL_QUAD
VISIT_CELL_TET = _simV2.VISIT_CELL_TET
VISIT_CELL_PYR = _simV2.VISIT_CELL_PYR
VISIT_CELL_WEDGE = _simV2.VISIT_CELL_WEDGE
VISIT_CELL_HEX = _simV2.VISIT_CELL_HEX
VISIT_CELL_POINT = _simV2.VISIT_CELL_POINT
VISIT_CELL_POLYHEDRON = _simV2.VISIT_CELL_POLYHEDRON
VISIT_COORD_MODE_SEPARATE = _simV2.VISIT_COORD_MODE_SEPARATE
VISIT_COORD_MODE_INTERLEAVED = _simV2.VISIT_COORD_MODE_INTERLEAVED
VISIT_GHOSTCELL_REAL = _simV2.VISIT_GHOSTCELL_REAL
VISIT_GHOSTCELL_INTERIOR_BOUNDARY = _simV2.VISIT_GHOSTCELL_INTERIOR_BOUNDARY
VISIT_GHOSTCELL_EXTERIOR_BOUNDARY = _simV2.VISIT_GHOSTCELL_EXTERIOR_BOUNDARY
VISIT_GHOSTCELL_ENHANCED_CONNECTIVITY = _simV2.VISIT_GHOSTCELL_ENHANCED_CONNECTIVITY
VISIT_GHOSTCELL_REDUCED_CONNECTIVITY = _simV2.VISIT_GHOSTCELL_REDUCED_CONNECTIVITY
VISIT_GHOSTCELL_BLANK = _simV2.VISIT_GHOSTCELL_BLANK
VISIT_CSG_QUADRIC_G = _simV2.VISIT_CSG_QUADRIC_G
VISIT_CSG_SPHERE_PR = _simV2.VISIT_CSG_SPHERE_PR
VISIT_CSG_ELLIPSOID_PRRR = _simV2.VISIT_CSG_ELLIPSOID_PRRR
VISIT_CSG_PLANE_G = _simV2.VISIT_CSG_PLANE_G
VISIT_CSG_PLANE_X = _simV2.VISIT_CSG_PLANE_X
VISIT_CSG_PLANE_Y = _simV2.VISIT_CSG_PLANE_Y
VISIT_CSG_PLANE_Z = _simV2.VISIT_CSG_PLANE_Z
VISIT_CSG_PLANE_PN = _simV2.VISIT_CSG_PLANE_PN
VISIT_CSG_PLANE_PPP = _simV2.VISIT_CSG_PLANE_PPP
VISIT_CSG_CYLINDER_PNLR = _simV2.VISIT_CSG_CYLINDER_PNLR
VISIT_CSG_CYLINDER_PPR = _simV2.VISIT_CSG_CYLINDER_PPR
VISIT_CSG_BOX_XYZXYZ = _simV2.VISIT_CSG_BOX_XYZXYZ
VISIT_CSG_CONE_PNLA = _simV2.VISIT_CSG_CONE_PNLA
VISIT_CSG_CONE_PPA = _simV2.VISIT_CSG_CONE_PPA
VISIT_CSG_POLYHEDRON_KF = _simV2.VISIT_CSG_POLYHEDRON_KF
VISIT_CSG_HEX_6F = _simV2.VISIT_CSG_HEX_6F
VISIT_CSG_TET_4F = _simV2.VISIT_CSG_TET_4F
VISIT_CSG_PYRAMID_5F = _simV2.VISIT_CSG_PYRAMID_5F
VISIT_CSG_PRISM_5F = _simV2.VISIT_CSG_PRISM_5F
VISIT_CSG_QUADRATIC_G = _simV2.VISIT_CSG_QUADRATIC_G
VISIT_CSG_CIRCLE_PR = _simV2.VISIT_CSG_CIRCLE_PR
VISIT_CSG_ELLIPSE_PRR = _simV2.VISIT_CSG_ELLIPSE_PRR
VISIT_CSG_LINE_G = _simV2.VISIT_CSG_LINE_G
VISIT_CSG_LINE_X = _simV2.VISIT_CSG_LINE_X
VISIT_CSG_LINE_Y = _simV2.VISIT_CSG_LINE_Y
VISIT_CSG_LINE_PN = _simV2.VISIT_CSG_LINE_PN
VISIT_CSG_LINE_PP = _simV2.VISIT_CSG_LINE_PP
VISIT_CSG_BOX_XYXY = _simV2.VISIT_CSG_BOX_XYXY
VISIT_CSG_ANGLE_PNLA = _simV2.VISIT_CSG_ANGLE_PNLA
VISIT_CSG_ANGLE_PPA = _simV2.VISIT_CSG_ANGLE_PPA
VISIT_CSG_POLYGON_KP = _simV2.VISIT_CSG_POLYGON_KP
VISIT_CSG_TRI_3P = _simV2.VISIT_CSG_TRI_3P
VISIT_CSG_QUAD_4P = _simV2.VISIT_CSG_QUAD_4P
VISIT_CSG_INNER = _simV2.VISIT_CSG_INNER
VISIT_CSG_OUTER = _simV2.VISIT_CSG_OUTER
VISIT_CSG_ON = _simV2.VISIT_CSG_ON
VISIT_CSG_UNION = _simV2.VISIT_CSG_UNION
VISIT_CSG_INTERSECT = _simV2.VISIT_CSG_INTERSECT
VISIT_CSG_DIFF = _simV2.VISIT_CSG_DIFF
VISIT_CSG_COMPLIMENT = _simV2.VISIT_CSG_COMPLIMENT
VISIT_CSG_XFORM = _simV2.VISIT_CSG_XFORM
VISIT_CSG_SWEEP = _simV2.VISIT_CSG_SWEEP
VISIT_IMAGEFORMAT_BMP = _simV2.VISIT_IMAGEFORMAT_BMP
VISIT_IMAGEFORMAT_JPEG = _simV2.VISIT_IMAGEFORMAT_JPEG
VISIT_IMAGEFORMAT_PNG = _simV2.VISIT_IMAGEFORMAT_PNG
VISIT_IMAGEFORMAT_POVRAY = _simV2.VISIT_IMAGEFORMAT_POVRAY
VISIT_IMAGEFORMAT_PPM = _simV2.VISIT_IMAGEFORMAT_PPM
VISIT_IMAGEFORMAT_RGB = _simV2.VISIT_IMAGEFORMAT_RGB
VISIT_IMAGEFORMAT_TIFF = _simV2.VISIT_IMAGEFORMAT_TIFF
VisIt_CommandMetaData_free = _simV2.VisIt_CommandMetaData_free
VisIt_CommandMetaData_setName = _simV2.VisIt_CommandMetaData_setName
VisIt_CSGMesh_free = _simV2.VisIt_CSGMesh_free
VisIt_CSGMesh_setRegions = _simV2.VisIt_CSGMesh_setRegions
VisIt_CSGMesh_setZonelist = _simV2.VisIt_CSGMesh_setZonelist
VisIt_CSGMesh_setBoundaryTypes = _simV2.VisIt_CSGMesh_setBoundaryTypes
VisIt_CSGMesh_setBoundaryCoeffs = _simV2.VisIt_CSGMesh_setBoundaryCoeffs
VisIt_CSGMesh_setExtents = _simV2.VisIt_CSGMesh_setExtents
VisIt_CurveData_free = _simV2.VisIt_CurveData_free
VisIt_CurveData_setCoordsXY = _simV2.VisIt_CurveData_setCoordsXY
VisIt_CurveMetaData_free = _simV2.VisIt_CurveMetaData_free
VisIt_CurveMetaData_setName = _simV2.VisIt_CurveMetaData_setName
VisIt_CurveMetaData_setXUnits = _simV2.VisIt_CurveMetaData_setXUnits
VisIt_CurveMetaData_setYUnits = _simV2.VisIt_CurveMetaData_setYUnits
VisIt_CurveMetaData_setXLabel = _simV2.VisIt_CurveMetaData_setXLabel
VisIt_CurveMetaData_setYLabel = _simV2.VisIt_CurveMetaData_setYLabel
VisIt_CurvilinearMesh_free = _simV2.VisIt_CurvilinearMesh_free
VisIt_CurvilinearMesh_setCoordsXY = _simV2.VisIt_CurvilinearMesh_setCoordsXY
VisIt_CurvilinearMesh_setCoordsXYZ = _simV2.VisIt_CurvilinearMesh_setCoordsXYZ
VisIt_CurvilinearMesh_setCoords2 = _simV2.VisIt_CurvilinearMesh_setCoords2
VisIt_CurvilinearMesh_setCoords3 = _simV2.VisIt_CurvilinearMesh_setCoords3
VisIt_CurvilinearMesh_setRealIndices = _simV2.VisIt_CurvilinearMesh_setRealIndices
VisIt_CurvilinearMesh_setBaseIndex = _simV2.VisIt_CurvilinearMesh_setBaseIndex
VisIt_CurvilinearMesh_setGhostCells = _simV2.VisIt_CurvilinearMesh_setGhostCells
VisIt_DomainBoundaries_free = _simV2.VisIt_DomainBoundaries_free
VisIt_DomainBoundaries_set_type = _simV2.VisIt_DomainBoundaries_set_type
VisIt_DomainBoundaries_set_numDomains = _simV2.VisIt_DomainBoundaries_set_numDomains
VisIt_DomainBoundaries_set_rectIndices = _simV2.VisIt_DomainBoundaries_set_rectIndices
VisIt_DomainBoundaries_set_amrIndices = _simV2.VisIt_DomainBoundaries_set_amrIndices
VisIt_DomainList_free = _simV2.VisIt_DomainList_free
VisIt_DomainList_setDomains = _simV2.VisIt_DomainList_setDomains
VisIt_DomainNesting_free = _simV2.VisIt_DomainNesting_free
VisIt_DomainNesting_set_dimensions = _simV2.VisIt_DomainNesting_set_dimensions
VisIt_DomainNesting_set_levelRefinement = _simV2.VisIt_DomainNesting_set_levelRefinement
VisIt_DomainNesting_set_nestingForPatch = _simV2.VisIt_DomainNesting_set_nestingForPatch
VisIt_ExpressionMetaData_free = _simV2.VisIt_ExpressionMetaData_free
VisIt_ExpressionMetaData_setName = _simV2.VisIt_ExpressionMetaData_setName
VisIt_ExpressionMetaData_setDefinition = _simV2.VisIt_ExpressionMetaData_setDefinition
VisIt_ExpressionMetaData_setType = _simV2.VisIt_ExpressionMetaData_setType
VisIt_MaterialData_free = _simV2.VisIt_MaterialData_free
VisIt_MaterialData_appendCells = _simV2.VisIt_MaterialData_appendCells
VisIt_MaterialData_addCleanCell = _simV2.VisIt_MaterialData_addCleanCell
VisIt_MaterialData_addMixedCell = _simV2.VisIt_MaterialData_addMixedCell
VisIt_MaterialData_setMaterials = _simV2.VisIt_MaterialData_setMaterials
VisIt_MaterialData_setMixedMaterials = _simV2.VisIt_MaterialData_setMixedMaterials
VisIt_MaterialMetaData_free = _simV2.VisIt_MaterialMetaData_free
VisIt_MaterialMetaData_setName = _simV2.VisIt_MaterialMetaData_setName
VisIt_MaterialMetaData_setMeshName = _simV2.VisIt_MaterialMetaData_setMeshName
VisIt_MaterialMetaData_addMaterialName = _simV2.VisIt_MaterialMetaData_addMaterialName
VisIt_MeshMetaData_free = _simV2.VisIt_MeshMetaData_free
VisIt_MeshMetaData_setName = _simV2.VisIt_MeshMetaData_setName
VisIt_MeshMetaData_setMeshType = _simV2.VisIt_MeshMetaData_setMeshType
VisIt_MeshMetaData_setTopologicalDimension = _simV2.VisIt_MeshMetaData_setTopologicalDimension
VisIt_MeshMetaData_setSpatialDimension = _simV2.VisIt_MeshMetaData_setSpatialDimension
VisIt_MeshMetaData_setNumDomains = _simV2.VisIt_MeshMetaData_setNumDomains
VisIt_MeshMetaData_setDomainTitle = _simV2.VisIt_MeshMetaData_setDomainTitle
VisIt_MeshMetaData_setDomainPieceName = _simV2.VisIt_MeshMetaData_setDomainPieceName
VisIt_MeshMetaData_addDomainName = _simV2.VisIt_MeshMetaData_addDomainName
VisIt_MeshMetaData_setNumGroups = _simV2.VisIt_MeshMetaData_setNumGroups
VisIt_MeshMetaData_setGroupTitle = _simV2.VisIt_MeshMetaData_setGroupTitle
VisIt_MeshMetaData_setGroupPieceName = _simV2.VisIt_MeshMetaData_setGroupPieceName
VisIt_MeshMetaData_addGroupId = _simV2.VisIt_MeshMetaData_addGroupId
VisIt_MeshMetaData_setXUnits = _simV2.VisIt_MeshMetaData_setXUnits
VisIt_MeshMetaData_setYUnits = _simV2.VisIt_MeshMetaData_setYUnits
VisIt_MeshMetaData_setZUnits = _simV2.VisIt_MeshMetaData_setZUnits
VisIt_MeshMetaData_setXLabel = _simV2.VisIt_MeshMetaData_setXLabel
VisIt_MeshMetaData_setYLabel = _simV2.VisIt_MeshMetaData_setYLabel
VisIt_MeshMetaData_setZLabel = _simV2.VisIt_MeshMetaData_setZLabel
VisIt_MeshMetaData_setCellOrigin = _simV2.VisIt_MeshMetaData_setCellOrigin
VisIt_MeshMetaData_setNodeOrigin = _simV2.VisIt_MeshMetaData_setNodeOrigin
VisIt_MeshMetaData_setHasSpatialExtents = _simV2.VisIt_MeshMetaData_setHasSpatialExtents
VisIt_MeshMetaData_setSpatialExtents = _simV2.VisIt_MeshMetaData_setSpatialExtents
VisIt_NameList_free = _simV2.VisIt_NameList_free
VisIt_NameList_addName = _simV2.VisIt_NameList_addName
VisIt_PointMesh_free = _simV2.VisIt_PointMesh_free
VisIt_PointMesh_setCoordsXY = _simV2.VisIt_PointMesh_setCoordsXY
VisIt_PointMesh_setCoordsXYZ = _simV2.VisIt_PointMesh_setCoordsXYZ
VisIt_PointMesh_setCoords = _simV2.VisIt_PointMesh_setCoords
VisIt_RectilinearMesh_free = _simV2.VisIt_RectilinearMesh_free
VisIt_RectilinearMesh_setCoordsXY = _simV2.VisIt_RectilinearMesh_setCoordsXY
VisIt_RectilinearMesh_setCoordsXYZ = _simV2.VisIt_RectilinearMesh_setCoordsXYZ
VisIt_RectilinearMesh_setRealIndices = _simV2.VisIt_RectilinearMesh_setRealIndices
VisIt_RectilinearMesh_setBaseIndex = _simV2.VisIt_RectilinearMesh_setBaseIndex
VisIt_RectilinearMesh_setGhostCells = _simV2.VisIt_RectilinearMesh_setGhostCells
VisIt_SimulationMetaData_free = _simV2.VisIt_SimulationMetaData_free
VisIt_SimulationMetaData_setMode = _simV2.VisIt_SimulationMetaData_setMode
VisIt_SimulationMetaData_setCycleTime = _simV2.VisIt_SimulationMetaData_setCycleTime
VisIt_SimulationMetaData_addMesh = _simV2.VisIt_SimulationMetaData_addMesh
VisIt_SimulationMetaData_addVariable = _simV2.VisIt_SimulationMetaData_addVariable
VisIt_SimulationMetaData_addMaterial = _simV2.VisIt_SimulationMetaData_addMaterial
VisIt_SimulationMetaData_addCurve = _simV2.VisIt_SimulationMetaData_addCurve
VisIt_SimulationMetaData_addExpression = _simV2.VisIt_SimulationMetaData_addExpression
VisIt_SimulationMetaData_addSpecies = _simV2.VisIt_SimulationMetaData_addSpecies
VisIt_SimulationMetaData_addGenericCommand = _simV2.VisIt_SimulationMetaData_addGenericCommand
VisIt_SimulationMetaData_addCustomCommand = _simV2.VisIt_SimulationMetaData_addCustomCommand
VisIt_SpeciesData_free = _simV2.VisIt_SpeciesData_free
VisIt_SpeciesData_addSpeciesName = _simV2.VisIt_SpeciesData_addSpeciesName
VisIt_SpeciesData_setSpecies = _simV2.VisIt_SpeciesData_setSpecies
VisIt_SpeciesData_setSpeciesMF = _simV2.VisIt_SpeciesData_setSpeciesMF
VisIt_SpeciesData_setMixedSpecies = _simV2.VisIt_SpeciesData_setMixedSpecies
VisIt_SpeciesMetaData_free = _simV2.VisIt_SpeciesMetaData_free
VisIt_SpeciesMetaData_setName = _simV2.VisIt_SpeciesMetaData_setName
VisIt_SpeciesMetaData_setMeshName = _simV2.VisIt_SpeciesMetaData_setMeshName
VisIt_SpeciesMetaData_setMaterialName = _simV2.VisIt_SpeciesMetaData_setMaterialName
VisIt_SpeciesMetaData_addSpeciesName = _simV2.VisIt_SpeciesMetaData_addSpeciesName
VisIt_UnstructuredMesh_free = _simV2.VisIt_UnstructuredMesh_free
VisIt_UnstructuredMesh_setCoordsXY = _simV2.VisIt_UnstructuredMesh_setCoordsXY
VisIt_UnstructuredMesh_setCoordsXYZ = _simV2.VisIt_UnstructuredMesh_setCoordsXYZ
VisIt_UnstructuredMesh_setCoords = _simV2.VisIt_UnstructuredMesh_setCoords
VisIt_UnstructuredMesh_setConnectivity = _simV2.VisIt_UnstructuredMesh_setConnectivity
VisIt_UnstructuredMesh_setRealIndices = _simV2.VisIt_UnstructuredMesh_setRealIndices
VisIt_UnstructuredMesh_setGhostCells = _simV2.VisIt_UnstructuredMesh_setGhostCells
VisIt_VariableData_free = _simV2.VisIt_VariableData_free
VisIt_VariableMetaData_free = _simV2.VisIt_VariableMetaData_free
VisIt_VariableMetaData_setName = _simV2.VisIt_VariableMetaData_setName
VisIt_VariableMetaData_setMeshName = _simV2.VisIt_VariableMetaData_setMeshName
VisIt_VariableMetaData_setUnits = _simV2.VisIt_VariableMetaData_setUnits
VisIt_VariableMetaData_setCentering = _simV2.VisIt_VariableMetaData_setCentering
VisIt_VariableMetaData_setType = _simV2.VisIt_VariableMetaData_setType
VisIt_VariableMetaData_setTreatAsASCII = _simV2.VisIt_VariableMetaData_setTreatAsASCII

VisItGetSockets = _simV2.VisItGetSockets
VisItReadConsole = _simV2.VisItReadConsole
VisIt_CommandMetaData_alloc = _simV2.VisIt_CommandMetaData_alloc
VisIt_CSGMesh_alloc = _simV2.VisIt_CSGMesh_alloc
VisIt_CurveData_alloc = _simV2.VisIt_CurveData_alloc
VisIt_CurveMetaData_alloc = _simV2.VisIt_CurveMetaData_alloc
VisIt_CurvilinearMesh_alloc = _simV2.VisIt_CurvilinearMesh_alloc
VisIt_DomainBoundaries_alloc = _simV2.VisIt_DomainBoundaries_alloc
VisIt_DomainList_alloc = _simV2.VisIt_DomainList_alloc
VisIt_DomainNesting_alloc = _simV2.VisIt_DomainNesting_alloc
VisIt_ExpressionMetaData_alloc = _simV2.VisIt_ExpressionMetaData_alloc
VisIt_MaterialData_alloc = _simV2.VisIt_MaterialData_alloc
VisIt_MaterialData_addMaterial = _simV2.VisIt_MaterialData_addMaterial
VisIt_MaterialMetaData_alloc = _simV2.VisIt_MaterialMetaData_alloc
VisIt_MeshMetaData_alloc = _simV2.VisIt_MeshMetaData_alloc
VisIt_NameList_alloc = _simV2.VisIt_NameList_alloc
VisIt_PointMesh_alloc = _simV2.VisIt_PointMesh_alloc
VisIt_RectilinearMesh_alloc = _simV2.VisIt_RectilinearMesh_alloc
VisIt_SimulationMetaData_alloc = _simV2.VisIt_SimulationMetaData_alloc
VisIt_SpeciesData_alloc = _simV2.VisIt_SpeciesData_alloc
VisIt_SpeciesMetaData_alloc = _simV2.VisIt_SpeciesMetaData_alloc
VisIt_UnstructuredMesh_alloc = _simV2.VisIt_UnstructuredMesh_alloc
VisIt_VariableData_alloc = _simV2.VisIt_VariableData_alloc
VisIt_VariableMetaData_alloc = _simV2.VisIt_VariableMetaData_alloc

