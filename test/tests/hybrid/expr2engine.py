# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  expr2engine.py
#
#  Tests:      mesh      - 2D, curvilinear, single domain
#              plots     - FilledBoundary
#              databases - PDB
#
#  Purpose:    This test case tests the viewer's ability to send not only the
#              user-defined expressions to the engine but also the correct
#              database expressions.
#
#  Programmer: Brad Whitlock
#  Date:       Fri Feb 18 14:01:48 PST 2005
#
#  Modifications:
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Cyrus Harrison, Thu Mar 25 09:57:34 PDT 2010
#    Added call(s) to DrawPlots() b/c of changes to the default plot state
#    behavior when an operator is added.
#
# ----------------------------------------------------------------------------

def TestExpressionList(name):
    exprList = Expressions()
    s = "Expressions:\n"
    index = 1
    for expr in exprList:
        s = s + "expression %d: %s = %s\n" % (index, expr[0], expr[1])
        index = index + 1
    TestText(name, s)

#
# Create some expressions.
#
DefineScalarExpression("user_defined1", "u * u")
DefineScalarExpression("user_defined2", "v + v")
DefineVectorExpression("user_defined3", "{u, v, w}")

# Open a database and make a plot.
OpenDatabase(silo_data_path("globe.silo"))

AddPlot("Vector", "vel")
v = VectorAttributes()
v.nVectors = 4000
SetPlotOptions(v)
DrawPlots()

v = GetView3D()
v.viewNormal = (-0.63515, 0.317784, 0.703987)
v.viewUp = (0.176786, 0.947058, -0.268008)
SetView3D(v)

# Test the image that we should have by this point. Also make sure that the
# expression list contains the database expressions for the first database.
Test("expr2engine_00")
TestExpressionList("expr2engine_01")

# Open a different database. The expression list should only contain the 
# database variables from the new database.
OpenDatabase(silo_data_path("noise.silo"))

TestExpressionList("expr2engine_02")

# Test that the plot from the old database, which was a plot of an expression
# from the first database can still be generated.
AddOperator("Transform")
DrawPlots()
Test("expr2engine_03")

Exit()
