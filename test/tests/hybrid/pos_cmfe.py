# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  pos_cmfe.py
#
#  Defect ID:  '5528, '5677
#
#  Programmer: Hank Childs
#  Date:       January 9, 2006
#
#  Modifications:
#    Jeremy Meredith, Tue Jul 15 10:43:58 EDT 2008
#    Changed number of vectors in vector plot to match the old behavior.
#    (We now account for how many domains there are.)
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Cyrus Harrison, Thu Mar 25 09:57:34 PDT 2010
#    Added call(s) to DrawPlots() b/c of changes to the default plot state
#    behavior when an operator is added.
#
# ----------------------------------------------------------------------------




OpenDatabase("../data/silo_%s_test_data/curv2d.silo"%SILO_MODE)


# Test general capability.
DefineScalarExpression("cmfe", "pos_cmfe(<../data/silo_%s_test_data/ucd2d.silo:d>, curvmesh2d, -1.)"%SILO_MODE)
AddPlot("Pseudocolor", "cmfe")
DrawPlots()
Test("pos_cmfe_01")

# Zonal vars
DefineScalarExpression("cmfeZ", "pos_cmfe(<../data/silo_%s_test_data/ucd2d.silo:p>, curvmesh2d, -1.)"%SILO_MODE)
DeleteAllPlots()
AddPlot("Pseudocolor", "cmfeZ")
DrawPlots()
Test("pos_cmfe_02")

# Test using the expression as a secondary variable.
ChangeActivePlotsVar("d")
AddOperator("Threshold")
t = ThresholdAttributes()
t.lowerBounds = (0.1)
t.listedVarNames = ("cmfe")
SetOperatorOptions(t)
DrawPlots()
Test("pos_cmfe_03")

# Test using cmfe in another expression.
DefineVectorExpression("disp", "{0, 3*cmfe}")
RemoveLastOperator()
d = DisplaceAttributes()
d.variable = "disp"
SetDefaultOperatorOptions(d)
AddOperator("Displace")
DrawPlots()
Test("pos_cmfe_04")

RemoveLastOperator()
ChangeActivePlotsVar("cmfe")
DrawPlots()

Query("MinMax")
t = GetQueryOutputString()
TestText("pos_cmfe_05", t)

# Now test error conditions.
DefineScalarExpression("cmfe2", "pos_cmfe(pressure, quadmesh)")
ChangeActivePlotsVar("cmfe2")
DrawPlots()
t = GetLastError()
TestText("pos_cmfe_06", t)

DefineScalarExpression("cmfe3", "pos_cmfe(<bad_file.silo:a_var>, curvmesh2d, 0.)")
ChangeActivePlotsVar("cmfe3")
DrawPlots()
t = GetLastError()
TestText("pos_cmfe_07", t)

DefineScalarExpression("cmfe4", "pos_cmfe(<../data/silo_%s_test_data/wave0100.silo:pressure>, curvmesh2d)"%SILO_MODE)
ChangeActivePlotsVar("cmfe4")
DrawPlots()
t = GetLastError()
TestText("pos_cmfe_08", t)

DefineScalarExpression("cmfe5", "pos_cmfe(<../data/silo_%s_test_data/wave0100.silo:pressure>)"%SILO_MODE)
ChangeActivePlotsVar("cmfe5")
DrawPlots()
t = GetLastError()
TestText("pos_cmfe_09", t)

# Invalid variable in new database.
DefineScalarExpression("cmfe6", "pos_cmfe(<../data/silo_%s_test_data/ucd2d.silo:xyz>, curvmesh2d, 0.)"%SILO_MODE)
ChangeActivePlotsVar("cmfe6")
DrawPlots()
t = GetLastError()
TestText("pos_cmfe_10", t)

# Enough error checking, let's do some more testing when it works.

# 3D, multi-block to multi-block.
DeleteAllPlots()
OpenDatabase("../data/silo_%s_test_data/multi_ucd3d.silo"%SILO_MODE)
DefineScalarExpression("cmfe7", "pos_cmfe(<../data/silo_%s_test_data/multi_curv3d.silo:d>, mesh1, 0.)"%SILO_MODE)
AddPlot("Pseudocolor", "cmfe7")
DrawPlots()
Test("pos_cmfe_11")

DeleteAllPlots()

# multi-block to single-block.
DeleteAllPlots()
OpenDatabase("../data/silo_%s_test_data/curv3d.silo"%SILO_MODE)
DefineScalarExpression("cmfe8", "pos_cmfe(<../data/silo_%s_test_data/multi_curv3d.silo:d>, curvmesh3d, 0.)"%SILO_MODE)
AddPlot("Pseudocolor", "cmfe8")
DrawPlots()
Test("pos_cmfe_12")

# single-block to multi-block.
DeleteAllPlots()
OpenDatabase("../data/silo_%s_test_data/multi_curv3d.silo"%SILO_MODE)
DefineScalarExpression("cmfe9", "pos_cmfe(<../data/silo_%s_test_data/curv3d.silo:d>, mesh1, 0.)"%SILO_MODE)
AddPlot("Pseudocolor", "cmfe9")
DrawPlots()
Test("pos_cmfe_13")

# single-block to multi-block -- vector
DeleteAllPlots()
OpenDatabase("../data/silo_%s_test_data/multi_ucd3d.silo"%SILO_MODE)
DefineVectorExpression("cmfe10", "pos_cmfe(gradient(<../data/silo_%s_test_data/curv3d.silo:d>), mesh1, {0,0,0})"%SILO_MODE)
AddPlot("Vector", "cmfe10")
vec = VectorAttributes()
vec.nVectors = 400*36
SetPlotOptions(vec)
DrawPlots()
Test("pos_cmfe_14")

Exit()
