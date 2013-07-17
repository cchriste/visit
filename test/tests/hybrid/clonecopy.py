# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  clonecopy.py
#
#  Tests:      the CloneWindow function.
#
#  Defect ID:  VisIt00003027
#
#  Programmer: Brad Whitlock
#  Date:       Wed Feb 12 11:59:14 PDT 2003
#
#  Modifications:
#    Kathleen Bonnell, Thu Aug 28 14:34:57 PDT 2003
#    Remove compound var name from subset plots.
#
#    Brad Whitlock, Wed Mar 9 09:15:30 PDT 2005
#    Removed deprecated functions.
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Kathleen Biagas, Thu Jul 11 08:16:36 PDT 2013
#    Removed legacy style annotation settings.
#
# ----------------------------------------------------------------------------

# Turn off annotation
a = AnnotationAttributes()
TurnOffAllAnnotations(a)
SetDefaultAnnotationAttributes(a)

# Set up a non-default annotation object.
a1 = AnnotationAttributes()
a1.axes2D.visible = 1
a1.axes2D.xAxis.label.visible = 0
a1.axes2D.yAxis.label.visible = 0
a1.axes2D.xAxis.title.visible = 0
a1.axes2D.yAxis.title.visible = 0
a1.axes3D.visible = 0
a1.axes3D.triadFlag = 0
a1.axes3D.bboxFlag = 0
a1.userInfoFlag = 0
a1.databaseInfoFlag = 0
a1.legendInfoFlag = 0
a1.foregroundColor = (0, 0, 0, 255)
a1.backgroundColor = (255, 255, 255, 255)
a1.gradientBackgroundStyle = a1.TopToBottom
a1.gradientColor1 = (204, 153, 255, 255)
a1.gradientColor2 = (0, 0, 0, 255)
a1.backgroundMode = a1.Gradient

# Set up the plots.
OpenDatabase(silo_data_path("wave.visit"))

AddPlot("Subset", "Material")
DrawPlots()
v = View3DAttributes()
v.viewNormal = (-0.427729, 0.776091, 0.463391)
v.focus = (4.37669, 0.376992, 2.57924)
v.viewUp = (0.67875, 0.614328, -0.402368)
v.viewAngle = 30.
v.parallelScale = 5.03337
v.nearPlane = -11.2758
v.farPlane = 11.2758
v.perspective = 1
SetView3D(v)
SetAnnotationAttributes(a1)

# Show what it looks like in window 1.
Test("clonecopy_Window1Setup1")
# Set the time state to later in time.
SetTimeSliderState(20)
Test("clonecopy_Window1Later")

# Create a few iterated cloned windows. This used to crash the viewer.
windowCount = 1
for i in range(3):
    CloneWindow()
    DrawPlots()
    Test("clonecopy_CloneOfWindow%d" % windowCount)
    windowCount = windowCount + 1

# Now that we're in another window that's the product of window cloning,
# make sure that we can set the active frame.
SetTimeSliderState(60)
Test("clonecopy_Window%dFrame60" % windowCount)

# Delete all but 2 windows
for i in range(3, windowCount+1):
    SetActiveWindow(i)
    DeleteWindow()

# Reset the windows that remain
SetActiveWindow(2)
ResetView()
DeleteAllPlots()
SetAnnotationAttributes(a)
SetActiveWindow(1)
ResetView()
DeleteAllPlots()
SetAnnotationAttributes(a)

# Add new plots
SetTimeSliderState(20)
AddPlot("Pseudocolor", "pressure")
AddPlot("Mesh", "quadmesh")
DrawPlots()
# Set the view
SetView3D(v)
# Make it have a gradient background color.
SetAnnotationAttributes(a1)
Test("clonecopy_Window1Setup2")

# Make window 2 the active window.
SetActiveWindow(2)

# Test CopyPlotsToWindow
CopyPlotsToWindow(1, 2)
DrawPlots()
Test("clonecopy_Window2CopyPlots")

# Test CopyViewToWindow
CopyViewToWindow(1, 2)
Test("clonecopy_Window2CopyView")

# Test CopyAnnotationsToWindow
CopyAnnotationsToWindow(1, 2)
Test("clonecopy_Window2CopyAnnotations")

Exit()
