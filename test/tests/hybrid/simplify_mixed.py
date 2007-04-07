# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  simplify_mixed.py
#
#  Tests:      plots     - filled boundary
#
#  Defect ID:  '4363, '6464, '6504
#
#  Programmer: Hank Childs
#  Date:       August 19, 2005
#
# ----------------------------------------------------------------------------

# Turn off all annotation
a = AnnotationAttributes()
a.axesFlag2D = 0
a.axesFlag = 0
a.triadFlag = 0
a.bboxFlag = 0
a.userInfoFlag = 0
a.databaseInfoFlag = 0
a.legendInfoFlag = 0
a.backgroundMode = 0
a.foregroundColor = (0, 0, 0, 255)
a.backgroundColor = (255, 255, 255, 255)
SetAnnotationAttributes(a)


OpenDatabase("../data/boxlib_test_data/2D/plt0822/Header")

AddPlot("FilledBoundary", "materials")
DrawPlots()

v = GetView2D()
v.windowCoords = (0.0084, 0.0215, 0.0920, 0.1034)
v.viewportCoords = (0.2, 0.95, 0.15, 0.95)
SetView2D(v)

Test("simplify_mixed01")

m = MaterialAttributes()
m.simplifyHeavilyMixedZones = 1
m.maxMaterialsPerZone = 2
SetMaterialAttributes(m)

AddPlot("Boundary", "materials")
b = BoundaryAttributes()
b.colorType = b.ColorBySingleColor
SetPlotOptions(b)

DrawPlots()

Test("simplify_mixed02")

DeleteAllPlots()
m.maxMaterialsPerZone = 1
SetMaterialAttributes(m)

AddPlot("FilledBoundary", "materials")
DrawPlots()

Test("simplify_mixed03")

Exit()
