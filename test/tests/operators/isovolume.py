# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  isovolume.py
#
#  Tests:      mesh      - 3D unstructured, single domain,
#                          3D rectilinear, multiple domain,
#                          2D curvilinear, multiple domain
#              plots     - pc, mesh, subset, contour
#              operators - isovolume
#              selection - none
#
#  Defect ID:  '5640
#
#  Programmer: Hank Childs
#  Date:       March 27, 2004
#
#  Modifications:
#
#    Hank Childs, Wed Nov 17 15:38:37 PST 2004
#    Added test for isovolumes of poly-data where the poly-data is totally
#    within the data range. ['5640]
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Cyrus Harrison, Thu Mar 25 09:57:34 PDT 2010
#    Added call(s) to DrawPlots() b/c of changes to the default plot state
#    behavior when an operator is added.
#
# ----------------------------------------------------------------------------


OpenDatabase(silo_data_path("globe.silo"))


AddPlot("Pseudocolor", "u")
DrawPlots()

v = GetView3D()
v.SetViewNormal(-0.528889, 0.367702, 0.7649)
v.SetViewUp(0.176641, 0.929226, -0.324558)
v.SetParallelScale(17.3205)
v.SetPerspective(1)
SetView3D(v)

#
# Normal PC isovolume for globe.
#
isovol = IsovolumeAttributes()

isovol.lbound = -4
isovol.ubound =  4
SetDefaultOperatorOptions(isovol)
AddOperator("Isovolume")
DrawPlots()
Test("ops_isovol01")


#
# Normal PC plot of globe, isovolumeing by a different variable than what we
# are coloring by.
#
RemoveAllOperators()
isovol.lbound = 140
isovol.ubound = 340
isovol.variable = "t"
SetDefaultOperatorOptions(isovol)
AddOperator("Isovolume")
DrawPlots()
Test("ops_isovol02")

DeleteAllPlots()

#
# Contour lines by one variable, isovolumeing by another.  Multi-block,
# curvilinear, 2D.
#
OpenDatabase(silo_data_path("multi_curv2d.silo"))

AddPlot("Contour", "u")
DrawPlots()

isovol.lbound = 0.7
isovol.ubound = 0.9
isovol.variable = "v"
SetDefaultOperatorOptions(isovol)
AddOperator("Isovolume")
DrawPlots()

Test("ops_isovol03")

DeleteAllPlots()

#
# Material plot, isovolumeed by a scalar variable.  Multi-block,
# curvilinear, 2D.
#
AddPlot("Subset", "mat1")
DrawPlots()

isovol.lbound = -0.4
isovol.ubound = 1.0
isovol.variable = "u"
SetDefaultOperatorOptions(isovol)
AddOperator("Isovolume")
DrawPlots()

Test("ops_isovol04")

DeleteAllPlots()

#
# Mesh plot and PC plot, both isovolumeed by the same variable, criteria.
# Multi-block, rectilinear, 3D.
#
OpenDatabase(silo_data_path("multi_rect3d.silo"))


isovol.lbound = 0.4
isovol.ubound = 1.0
isovol.variable = "u"
SetDefaultOperatorOptions(isovol)

AddPlot("Pseudocolor", "u")
AddOperator("Isovolume")
AddPlot("Mesh", "mesh1")
AddOperator("Isovolume")
DrawPlots()

Test("ops_isovol05")

#
# Test that we can slice poly-data correctly.
#
DeleteAllPlots()
OpenDatabase(silo_data_path("rect3d.silo"))


isovol.lbound = -1.0
isovol.ubound = 1.0
isovol.variable = "default"
SetDefaultOperatorOptions(isovol)

AddPlot("Pseudocolor", "d")
AddOperator("Slice")
AddOperator("Isovolume")
DrawPlots()

Test("ops_isovol06")

isovol.lbound = 0.5
SetOperatorOptions(isovol)
Test("ops_isovol07")

Exit()
