# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  expr_cmfe.py
#
#  Defect ID:  None
#
#  Programmer: Hank Childs
#  Date:       September 9, 2005
#
#  Modifications:
#
#    Hank Childs, Thu Dec 29 11:21:26 PST 2005
#    Expand the color range for a plot of an expression that should result in 
#    uniformly "1", but actually has some small variation.  This causes 
#    issues with coloring between optimized and non-optimized modes.
#
# ----------------------------------------------------------------------------



OpenDatabase("../data/silo_hdf5_test_data/wave.visit")


# Test that database expressions can still be generated.
DefineVectorExpression("cmfe", "conn_cmfe(<../data/silo_hdf5_test_data/wave.visit[30]i:direction>, quadmesh)")
AddPlot("Vector", "cmfe")
DrawPlots()
Test("expr_cmfe_01")

DeleteAllPlots()
DefineScalarExpression("cmfe2", "conn_cmfe(coord(<../data/silo_hdf5_test_data/wave.visit[40]i:pressure>)[1], quadmesh)")
AddPlot("Pseudocolor", "cmfe2")
DrawPlots()
Test("expr_cmfe_02")

DeleteAllPlots()
DefineScalarExpression("cmfe3", "coord(quadmesh)[1] - conn_cmfe(coord(<../data/silo_hdf5_test_data/wave.visit[40]i:pressure>)[1], quadmesh)")
AddPlot("Pseudocolor", "cmfe3")
DrawPlots()
Test("expr_cmfe_03")

DeleteAllPlots()
DefineScalarExpression("cmfe4", "coord(quadmesh)[1] - cmfe2")
AddPlot("Pseudocolor", "cmfe4")
DrawPlots()
Test("expr_cmfe_04")

DeleteAllPlots()
DefineScalarExpression("cmfe5", "volume(quadmesh) / conn_cmfe(volume(<../data/silo_hdf5_test_data/wave.visit[40]i:quadmesh>), quadmesh)")
AddPlot("Pseudocolor", "cmfe5")
pc = PseudocolorAttributes()
pc.min = 0.5
pc.minFlag = 1
pc.max = 1.5
pc.maxFlag = 1
SetPlotOptions(pc)
DrawPlots()
Test("expr_cmfe_05")


Exit()
