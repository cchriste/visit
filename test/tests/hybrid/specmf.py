# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  val4mat.py
#
#  Tests:      mesh      - 2d structured
#              plots     - pc
#
#  Notes
#
#  Programmer: Cyrus Harrison
#  Date:       Tuesday 12, 2008
#
#  Modifiations:
#
# ----------------------------------------------------------------------------


OpenDatabase("../data/silo_%s_test_data/specmix_quad.silo"%SILO_MODE)
atts = PseudocolorAttributes()
atts.minFlag = 1
atts.min = 0.0
atts.maxFlag = 1
atts.max = 1.0
SetDefaultPlotOptions(atts)

# view the per material values for each of the 3 materials

DefineScalarExpression("spec_mix", "specmf(Species,1,1)")
AddPlot("Pseudocolor", "spec_mix")
DrawPlots()
Test("specmf_0")

OpenDatabase("../data/silo_%s_test_data/specmix_double_quad.silo"%SILO_MODE)
AddPlot("Pseudocolor", "spec_mix")
DrawPlots()
Test("specmf_1")

Exit()
