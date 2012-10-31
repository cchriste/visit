# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  box.py
#
#  Tests:      mesh      - 3D point, single domain,
#                          3D rectilinear, multiple domain.
#                          3D unstructured, single domain.
#              plots     - pc, vector, contour, vector
#
#  Defect ID:  3160
#
#  Programmer: Hank Childs
#  Date:       March 26, 2003
#
#  Modifications:
#    Kathleen Bonnell, Thu Aug 28 14:34:57 PDT 2003
#    Remove compound var name from subset plots.
#
#    Hank Childs, Wed Jun  4 08:56:08 PDT 2008
#    Test facelist filter after applying a box.
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
# ----------------------------------------------------------------------------


# 3D, rectilinear.  Multi-block  Tests ghost zones as well.
OpenDatabase(silo_data_path("bigsil.silo"))


atts = BoxAttributes()
atts.amount = 0
atts.minx = 0.4
atts.maxx = 0.7
atts.miny = 0.1
atts.maxy = 0.9
atts.minz = 0.4
atts.maxz = 1.7

AddPlot("Pseudocolor", "dist")
AddOperator("Box")
SetOperatorOptions(atts)
DrawPlots()
view = GetView3D()
view.viewNormal = (0.79, -0.035, 0.61)
view.focus = (0.5, 0.5, 0.5)
view.viewUp = (0, 1, 0)
view.parallelScale = 0.765
view.nearPlane = -0.25
view.farPlane = 1.73
SetView3D(view)
Test("ops_box01")
DeleteAllPlots()


# 3D, unstructured
OpenDatabase(silo_data_path("globe.silo"))

AddPlot("Subset", "mat1")
AddOperator("Box")
atts.minx = -5
atts.maxx = +5
atts.miny = -5
atts.maxy = +5
atts.minz = -5
atts.maxz = +5
SetOperatorOptions(atts)
DrawPlots()
view.viewNormal = (0.56, 0.81, -0.14)
view.focus = (0, 0, 0)
view.viewUp = (-0.66, 0.55, 0.49)
view.parallelScale = 11
view.nearPlane = -34
view.farPlane = 34
SetView3D(view)
Test("ops_box02")


atts.amount = 1
SetOperatorOptions(atts)
Test("ops_box03")

DeleteAllPlots()

view.parallelScale = 14
SetView3D(view)

OpenDatabase(silo_data_path("noise.silo"))

AddPlot("Vector", "airVfGradient")
vec = VectorAttributes()
SetPlotOptions(vec)
atts.minx = -5
atts.maxx = +5
atts.miny = -10
atts.maxy = +10
atts.minz = -10
atts.maxz = +10
AddOperator("Box")
SetOperatorOptions(atts)
DrawPlots()
Test("ops_box04")

DeleteAllPlots()

AddPlot("Contour", "hardyglobal")
AddOperator("Box")
SetOperatorOptions(atts)
DrawPlots()
Test("ops_box05")

DeleteAllPlots()

AddPlot("Mesh", "PointMesh")
mesh = MeshAttributes()
mesh.pointType = mesh.Box
mesh.pointSizeVar = "PointVar"
mesh.pointSizeVarEnabled = 1
SetPlotOptions(mesh)
AddOperator("Box")
SetOperatorOptions(atts)
DrawPlots()
Test("ops_box06")

# This final test should be blank.  The box operator was crashing with
# point meshes that fell outside the box previously.  ('3160)
atts.minx = -30
atts.maxx = -20
SetOperatorOptions(atts)
Test("ops_box07")

# If we apply the box operator to a rectilinear grid, it can change
# the "Extents", which can screw up the facelist filter.  Test that
# no one has undone the fix from Paul Selby.
DeleteAllPlots()
OpenDatabase(silo_data_path("rect3d.silo"))

AddPlot("Pseudocolor", "d")
AddOperator("Box")
atts.minx = 0.2
atts.maxx = 1.0
SetOperatorOptions(atts)
DrawPlots()
view.focus = (0.5, 0.5, 0.5)
view.parallelScale = 1.4
SetView3D(view)
Test("ops_box08")

Exit()
