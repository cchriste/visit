# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  singlecolor.py
#
#  Tests:      mesh      - 3D rectilinear, single domain
#              plots     - Boundary, FilledBoundary, Subset
#
#  Defect ID:  VisIt00002372
#
#  Programmer: Brad Whitlock
#  Date:       Thu Oct 23 15:36:30 PST 2003
#
#  Modifications:
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
# ----------------------------------------------------------------------------

# Set the single color to light blue and partially transparent using the
# plot's global opacity setting.
def SetSingleColor(atts):
    atts.singleColor = (153, 204, 255, 255)
    atts.colorType = b.ColorBySingleColor
    atts.opacity = 0.4
    SetPlotOptions(atts)

# Set the view that we want to use.
def InitializeView():
    v = View3DAttributes()
    v.viewNormal = (-0.428395 ,0.549517, 0.717293)
    v.focus = (0.5, 0.5, 0.5)
    v.viewUp = (0.186332, 0.830487, -0.52495)
    v.viewAngle = 30
    v.parallelScale = 0.866025
    v.nearPlane = -1.73205
    v.farPlane = 1.73205
    v.imagePan = (0.0183269, -0.0257188)
    v.imageZoom = 1.17591
    v.perspective = 1
    v.eyeAngle = 2
    SetView3D(v)

# Open the database.
OpenDatabase(silo_data_path("rect3d.silo"))


# Test the single color opacity for the Boundary plot
AddPlot("Boundary", "mat1")
b = BoundaryAttributes()
SetSingleColor(b)
DrawPlots()
InitializeView()
Test("singlecolor00")

# Test the single color opacity for the FilledBoundary plot
DeleteAllPlots()
AddPlot("FilledBoundary", "mat1")
f = FilledBoundaryAttributes()
SetSingleColor(f)
DrawPlots()
Test("singlecolor01")

# Test the single color opacity for the Subset plot
DeleteAllPlots()
AddPlot("Subset", "mat1")
s = SubsetAttributes()
SetSingleColor(s)
DrawPlots()
Test("singlecolor02")

Exit()
