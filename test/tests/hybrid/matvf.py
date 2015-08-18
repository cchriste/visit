# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  expressions.py
#
#  Tests:      mesh      - 3D unstructured, multi- domain
#                          3D rectilinear, single domain
#              plots     - Pseudocolor
#              operators - Threshold
#              selection - Material
#
#  Defect ID:  '2415 '3939 '5521
#
#  Programmer: Hank Childs
#  Date:       August 13, 2003
#
#  Modifications:
#
#    Hank Childs, Mon Oct 20 16:32:46 PDT 2003
#    Added bigsil portion to test '3939 (matvf with databases that create
#    ghost zones).
#
#    Hank Childs, Mon Jan  3 09:28:39 PST 2005
#    Added tests for materror.
#
#    Cyrus Harrison, Tue Feb 12 15:16:36 PST 2008
#    Added explicit tests for matvf w/ ghost zones. 
#
#    Cyrus Harrison, Wed Apr  9 11:17:22 PDT 2008
#    Added new test for post ghost info.
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Cyrus Harrison, Thu Mar 25 09:57:34 PDT 2010
#    Added call(s) to DrawPlots() b/c of changes to the default plot state
#    behavior when an operator is added.
#
# ----------------------------------------------------------------------------



OpenDatabase(silo_data_path("rect3d.silo"))



DefineScalarExpression("m1", 'matvf(mat1, ["1", "3"])')
AddPlot("Pseudocolor", "m1")
DrawPlots()

v = View3DAttributes()
v.viewNormal = (0.17, 0.74, 0.64)
v.focus = (0.5, 0.5, 0.5)
v.viewUp = (-0.31, 0.65, -0.68)
v.parallelScale = 0.866
v.nearPlane = -1.73
v.farPlane = 1.73
SetView3D(v)

Test("matvf_01")

t = ThresholdAttributes()
t.lowerBounds = (0.01)
t.upperBounds = (0.99)
SetDefaultOperatorOptions(t)
AddOperator("Threshold")
DrawPlots()

Test("matvf_02")

sil = SILRestriction()
matsets = sil.SetsInCategory("mat1")
sil.TurnOffAll()
sil.TurnOnSet(matsets[1])
SetPlotSILRestriction(sil)

Test("matvf_03")

DeleteAllPlots()


OpenDatabase(silo_data_path("multi_ucd3d.silo"))

AddPlot("Pseudocolor", "m1")
AddOperator("Threshold")
DrawPlots()

v = View3DAttributes()
v.viewNormal = (-0.82, 0.44, 0.37)
v.focus = (0, 3.3, 10)
v.viewUp = (0.13, 0.76, -0.62)
v.parallelScale = 7.32
v.nearPlane = -14.65
v.farPlane = 14.65
SetView3D(v)

Test("matvf_04")

sil = SILRestriction()
matsets = sil.SetsInCategory("mat1")
sil.TurnOffAll()
sil.TurnOnSet(matsets[1])
SetPlotSILRestriction(sil)

Test("matvf_05")

DeleteAllPlots()

DefineScalarExpression("mbs", 'matvf(mat, 1)')
OpenDatabase(silo_data_path("bigsil.silo"))

AddPlot("Pseudocolor", "mbs")
t = ThresholdAttributes()
t.lowerBounds = (0.5)
t.upperBounds = (1.1)
SetDefaultOperatorOptions(t)
AddOperator("Threshold")
DrawPlots()

v = View3DAttributes()
v.viewNormal = (-0.59, 0.33, 0.74)
v.focus = (0.5, 0.5, 0.5)
v.viewUp = (0.18, 0.94, -0.28)
v.parallelScale = 0.866
v.nearPlane = -1.73
v.farPlane = 1.73
SetView3D(v)

v = View3DAttributes()
v.viewNormal = (0.17, 0.74, 0.64)
v.focus = (0.5, 0.5, 0.5)
v.viewUp = (-0.31, 0.65, -0.68)
v.parallelScale = 0.866
v.nearPlane = -1.73
v.farPlane = 1.73
SetView3D(v)

Test("matvf_06")

#NOTE: the fix for '3939 included disabling the communication of ghost zones. 
# As a result, the matvf_07 picture will contain internal surfaces that should
# not be there.  If this ever gets fixed, then the baseline picture will need
# to be reset.
v.nearPlane = -0.1
SetView3D(v)

Test("matvf_07")

# Now test the materror expression.
DeleteAllPlots()
DefineScalarExpression("materror1", 'materror(mat1, 1)')
OpenDatabase(silo_data_path("rect3d.silo"))

AddPlot("Pseudocolor", "materror1")
sil = SILRestriction()
sil.TurnOnAll()
SetPlotSILRestriction(sil)
DrawPlots()

v = View3DAttributes()
v.viewNormal = (0.17, 0.74, 0.64)
v.focus = (0.5, 0.5, 0.5)
v.viewUp = (-0.31, 0.65, -0.68)
v.parallelScale = 0.866
v.nearPlane = -1.73
v.farPlane = 1.73
SetView3D(v)

Test("matvf_08")

DefineScalarExpression("abs_materror1", 'abs(materror1)')
t = ThresholdAttributes()
t.lowerBounds = (0.01)
t.upperBounds = (1.01)
t.listedVarNames = ("abs_materror1")
SetDefaultOperatorOptions(t)
AddOperator("Threshold")
DrawPlots()

Test("matvf_09")

# Now test that it can calculate "volume fractions" using areas.
ReplaceDatabase(silo_data_path("rect2d.silo"))

RemoveLastOperator()
Test("matvf_10")


DeleteAllPlots()
CloseDatabase(silo_data_path("rect2d.silo"))

OpenDatabase(silo_data_path("thinplane.silo"))


# tests for matvf w/ ghost zones

# test standard
DefineScalarExpression("mvf_2", "matvf(mat,2)")
AddPlot("Pseudocolor", "mvf_2")
atts = PseudocolorAttributes(1)
atts.minFlag = 1
atts.min = 0.0
atts.maxFlag = 1
atts.max = 1.0
atts.centering = atts.Natural
SetPlotOptions(atts)
DrawPlots()
Test("matvf_post_ghost_1")

# test with ghosts
atts = PseudocolorAttributes(1)
atts.minFlag = 1
atts.min = 0.0
atts.maxFlag = 1
atts.max = 1.0
atts.centering = atts.Nodal
SetPlotOptions(atts)
DrawPlots()
Test("matvf_post_ghost_2")

# switch back to ensure w/ ghost case caching does not alter w/o ghosts
atts = PseudocolorAttributes(1)
atts.minFlag = 1
atts.min = 0.0
atts.maxFlag = 1
atts.max = 1.0
atts.centering = atts.Zonal
SetPlotOptions(atts)
DrawPlots()
Test("matvf_post_ghost_3")

DeleteAllPlots()
CloseDatabase(silo_data_path("thinplane.silo"))

OpenDatabase(silo_data_path("bigsil.silo"))


# test for gradient w/ matvf, to check if post ghosts were correctly 
# requested. 

DefineScalarExpression("gate", "agrad(dist)[0] * matvf(mat,1)")
AddPlot("Pseudocolor", "gate")
DrawPlots()
Test("matvf_post_ghost_4")
DeleteAllPlots()


Exit()
