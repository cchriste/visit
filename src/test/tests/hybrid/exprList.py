# ---------------------------------------------------------------------------- 
#  CLASSES: nightly
#
#  Test Case:  exprList.py
#  Tests:      Expression list contents as windows are added and we switch
#              between databases that have expressions.
#
#  Defect ID:  VisIt00003955
#
#  Programmer: Brad Whitlock
#  Date:       Fri Oct 24 18:06:01 PST 2003
#
#  Modifications:
#    
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
# ----------------------------------------------------------------------------

TurnOnAllAnnotations()

# Define some expressions just so we have some in the list.
DefineScalarExpression("var1", "var2 + var3")
DefineScalarExpression("var4", "var5 * var6")
DefineScalarExpression("var7", "var8 / var9")

# Open the first database, which has some expressions.
OpenDatabase(silo_data_path("globe.silo"))

AddPlot("Pseudocolor", "speed")
DrawPlots()

# This test should show our scalar expressions + globe's expressions.
TestExpressions("exprList00")

# Add a new window and open a different database that has no expressions of
# its own.
AddWindow()
SetActiveWindow(2)
DeleteAllPlots()
OpenDatabase(silo_data_path("wave*.silo database"))

TestExpressions("exprList01")

# Going back to window 1, where globe is open. This should make the expression
# list contain globe's expressions.
SetActiveWindow(1)
TestExpressions("exprList02")

# Open a new database. This should make the expression list contain the
# expressions for rect3d and our scalar expressions.
OpenDatabase(silo_data_path("rect3d.silo"))

TestExpressions("exprList03")

# Add a plot
AddPlot("Pseudocolor", "u")
DrawPlots()

# Make the active plot be the plot of globe. The expression list should 
# contain globe's expressions.
SetActivePlots(0)
TestExpressions("exprList04")

Exit()
