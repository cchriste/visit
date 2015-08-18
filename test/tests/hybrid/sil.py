# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  sil.py
#
#  Tests:      mesh      - 3D curvilinear,multi-domain,ghost zones replicated.
#              plots     - mat subset, domain subset
#
#  Defect ID:  none
#
#  Programmer: Hank Childs
#  Date:       December 5, 2002
#
#  Modifications:
#    Kathleen Bonnell, Thu Aug 28 14:34:57 PDT 2003
#    Remove compound var name from subset plots.
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Hank Childs, Fri Feb 24 15:45:41 PST 2012
#    Add test for preserving SILs across time slice changes.
#
# ----------------------------------------------------------------------------


view = View3DAttributes()
view.viewNormal = (0.557976, 0.651128, 0.514485)
view.focus = (0.5, 0.5, 0.5)
view.viewUp = (-0.0955897, 0.666272, -0.739557)
view.viewAngle = 30
view.parallelScale = 0.866025
view.nearPlane = -1.73205
view.farPlane = 1.73205
view.perspective = 1
SetView3D(view)

OpenDatabase(silo_data_path("bigsil.silo"))

AddPlot("Subset", "mat")
DrawPlots()

# Test the normal material plot.
Test("sil1")

# Make sure that the ghost zones were generated correctly.
view.nearPlane = -0.3
SetView3D(view)
Test("sil2")

view.nearPlane = -1.73205
SetView3D(view)

TurnMaterialsOff("1")
Test("sil3")

TurnMaterialsOff()
TurnMaterialsOn("1")
Test("sil4")

DeleteAllPlots()

# Test that the SIL from the previous plot is preserved.
AddPlot("Subset", "domains")
DrawPlots()
Test("sil5")

OpenDatabase(data_path("Chombo_test_data/chombo.visit"))

AddPlot("Pseudocolor", "Density")
s = SILRestriction()
s.TurnOffSet(s.SetsInCategory("materials")[1])
SetPlotSILRestriction(s)
DrawPlots()
TimeSliderSetState(4)
s = SILRestriction()
if (s.UsesData(s.SetsInCategory("materials")[1])):
   str="Material 1 got turned back on!  (incorrect)\n"
else:
   str="Material 1 was correctly left off.\n"

TestText("sil6", str)
Exit()


