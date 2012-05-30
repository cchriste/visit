# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  ddf.py
#
#  Defect ID:  '5203
#
#  Programmer: Hank Childs
#  Date:       February 20, 2006
#
#  Modifications:
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
# ----------------------------------------------------------------------------



OpenDatabase(silo_data_path("curv2d.silo"))



AddPlot("Pseudocolor", "d")
DrawPlots()

t = ConstructDDFAttributes()
t.ddfName = "ddf1"
t.varnames = ("u")
t.ranges = (-1, 1)
t.numSamples = (4)
t.codomainName = "u"
t.statisticalOperator = t.Average
ConstructDDF(t)

DefineScalarExpression("e1", "u - apply_ddf(curvmesh2d, ddf1)")
ChangeActivePlotsVar("e1")
Test("ddf_01")

t.ddfName = "ddf2"
t.statisticalOperator = t.Maximum
t.codomainName = "v"
t.varnames = ("v")
ConstructDDF(t)
DefineScalarExpression("e2", "v - apply_ddf(curvmesh2d, ddf2)")
ChangeActivePlotsVar("e2")
Test("ddf_02")

t.ddfName = "ddf3"
t.varnames = ("u", "v")
t.ranges = (-1, 1, -1, 1)
t.numSamples = (25, 25)
t.codomainName = "u"
t.statisticalOperator = t.Minimum
ConstructDDF(t)

DefineScalarExpression("e3", "u - apply_ddf(curvmesh2d, ddf3)")
ChangeActivePlotsVar("e3")
Test("ddf_03")


ChangeActivePlotsVar("u")
t.ddfName = "ddf4"
t.varnames = ("u", "v")
t.ranges = (-1, 1, -1, 1)
t.numSamples = (25, 25)
t.codomainName = "u"
t.statisticalOperator = t.RMS
ConstructDDF(t)

DefineScalarExpression("e4", "apply_ddf(curvmesh2d, ddf4)")
ChangeActivePlotsVar("e4")
Test("ddf_04")

Exit()
