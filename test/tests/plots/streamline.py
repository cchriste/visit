# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  streamline.py
#
#  Tests:      streamline - 3D rectilinear
#              plots      - streamline
#
#  Defect ID:  -
#
#  Programmer: Christoph Garth
#  Date:       July 20, 2010
#
# ----------------------------------------------------------------------------

OpenDatabase("../data/silo_%s_test_data/noise.silo"%SILO_MODE)

AddPlot("Streamline", "grad", 1, 0)

View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0.270729, 0.624198, 0.732859)
View3DAtts.focus = (0.496062, 0.99603, 0.496062)
View3DAtts.viewUp = (-0.0922782, 0.774611, -0.62567)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 12.1829
View3DAtts.nearPlane = -24.3658
View3DAtts.farPlane = 24.3658
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
SetView3D(View3DAtts)

StreamlineAtts = StreamlineAttributes()
StreamlineAtts.sourceType = StreamlineAtts.SpecifiedPlane  # SpecifiedPoint, SpecifiedPointList, SpecifiedLine, SpecifiedCircle, SpecifiedPlane, SpecifiedSphere, SpecifiedBox
StreamlineAtts.maxStepLength = 0.02
StreamlineAtts.termination = 200
StreamlineAtts.pointSource = (0, 0, 0)
StreamlineAtts.lineStart = (0, 0, 0)
StreamlineAtts.lineEnd = (1, 0, 0)
StreamlineAtts.planeOrigin = (0.5, 1, 0.5)
StreamlineAtts.planeNormal = (0, 1, 0)
StreamlineAtts.planeUpAxis = (1, 0, 0)
StreamlineAtts.radius = 1
StreamlineAtts.sphereOrigin = (0, 0, 0)
StreamlineAtts.boxExtents = (0, 1, 0, 1, 0, 1)
StreamlineAtts.useWholeBox = 1
StreamlineAtts.pointList = (0, 0, 0, 1, 0, 0, 0, 1, 0)
StreamlineAtts.sampleDensity0 = 5
StreamlineAtts.sampleDensity1 = 5
StreamlineAtts.sampleDensity2 = 1
StreamlineAtts.displayMethod = StreamlineAtts.Tubes  # Lines, Tubes, Ribbons
StreamlineAtts.showSeeds = 1
StreamlineAtts.showHeads = 0
StreamlineAtts.tubeRadius = 0.25
StreamlineAtts.ribbonWidth = 0.125
StreamlineAtts.lineWidth = 2
StreamlineAtts.coloringMethod = StreamlineAtts.ColorBySpeed  # Solid, ColorBySpeed, ColorByVorticity, ColorByLength, ColorByTime, ColorBySeedPointID, ColorByVariable
StreamlineAtts.colorTableName = "Default"
StreamlineAtts.singleColor = (255, 0, 0, 255)
StreamlineAtts.legendFlag = 1
StreamlineAtts.lightingFlag = 1
StreamlineAtts.streamlineDirection = StreamlineAtts.Both  # Forward, Backward, Both
StreamlineAtts.relTol = 1e-06
StreamlineAtts.absTol = 1e-07
StreamlineAtts.terminationType = StreamlineAtts.Time  # Distance, Time, Step
StreamlineAtts.integrationType = StreamlineAtts.DormandPrince  # DormandPrince, AdamsBashforth, M3DC1Integrator
StreamlineAtts.streamlineAlgorithmType = StreamlineAtts.ParallelStaticDomains  # LoadOnDemand, ParallelStaticDomains, MasterSlave
StreamlineAtts.maxStreamlineProcessCount = 10
StreamlineAtts.maxDomainCacheSize = 3
StreamlineAtts.workGroupSize = 32
StreamlineAtts.pathlines = 0
StreamlineAtts.coloringVariable = "hardyglobal"
StreamlineAtts.legendMinFlag = 0
StreamlineAtts.legendMaxFlag = 0
StreamlineAtts.legendMin = 0
StreamlineAtts.legendMax = 1
StreamlineAtts.displayBegin = 0
StreamlineAtts.displayEnd = 1
StreamlineAtts.displayBeginFlag = 0
StreamlineAtts.displayEndFlag = 0
StreamlineAtts.seedDisplayRadius = 0.8
StreamlineAtts.headDisplayType = StreamlineAtts.Sphere  # Sphere, Cone
StreamlineAtts.headDisplayRadius = 0.25
StreamlineAtts.headDisplayHeight = 0.5
StreamlineAtts.opacityType = StreamlineAtts.None  # None, Constant, Ramp, VariableRange
StreamlineAtts.opacityVariable = ""
StreamlineAtts.opacity = 1
StreamlineAtts.opacityVarMin = 0
StreamlineAtts.opacityVarMax = 1
StreamlineAtts.opacityVarMinFlag = 0
StreamlineAtts.opacityVarMaxFlag = 0
StreamlineAtts.tubeDisplayDensity = 10
StreamlineAtts.geomDisplayQuality = StreamlineAtts.Medium  # Low, Medium, High, Super
StreamlineAtts.sampleDistance0 = 18
StreamlineAtts.sampleDistance1 = 18
StreamlineAtts.sampleDistance2 = 1
StreamlineAtts.fillInterior = 1
StreamlineAtts.randomSamples = 0
StreamlineAtts.randomSeed = 0
StreamlineAtts.numberOfRandomSamples = 1
StreamlineAtts.forceNodeCenteredData = 0

# test coloring options

StreamlineAtts.coloringMethod = StreamlineAtts.ColorBySpeed
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_01" )

StreamlineAtts.coloringMethod = StreamlineAtts.ColorByVorticity
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_02" )

StreamlineAtts.coloringMethod = StreamlineAtts.ColorByLength
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_03" )

StreamlineAtts.coloringMethod = StreamlineAtts.ColorByTime
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_04" )

StreamlineAtts.coloringMethod = StreamlineAtts.ColorBySeedPointID
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_05" )

StreamlineAtts.coloringMethod = StreamlineAtts.ColorByVariable
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_06" )

StreamlineAtts.coloringMethod = StreamlineAtts.Solid
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_07" )

# test termination modes (termination by time implied in previous tests)

StreamlineAtts.coloringMethod = StreamlineAtts.Solid

StreamlineAtts.terminationType = StreamlineAtts.Distance
StreamlineAtts.termination = 10
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_08" )

StreamlineAtts.termination = 5
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_09" )

StreamlineAtts.terminationType = StreamlineAtts.Step
StreamlineAtts.termination = 500
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_10" )

StreamlineAtts.terminationType = StreamlineAtts.Step
StreamlineAtts.termination = 250
SetPlotOptions(StreamlineAtts)
DrawPlots()
Test( "streamlines_11" )

Exit()
