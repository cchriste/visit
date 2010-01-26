# ----------------------------------------------------------------------------
#  MODES: serial
#  CLASSES: nightly
#
#  Test Case:  axistitles.py 
#
#  Tests:      Tests setting axis titles and units.
#
#  Programmer: Brad Whitlock
#  Date:       Thu Jul 28 11:07:57 PDT 2005
#
#  Modifications:
#    Brad Whitlock, Wed Apr 2 16:44:26 PST 2008
#    Modified the 3D test since setting the 3D font scale now actually works.
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
# ----------------------------------------------------------------------------

def SaveTestImage(name):
    # Save these images somewhat larger than a regular test case image
    # since the images contain a lot of text.
    swa = SaveWindowAttributes()
    swa.width = 500
    swa.height = 500
    swa.screenCapture = 0
    Test(name, swa)

#
# Test replacing 2D titles and units.
#
def Test2D():
    TestSection("Setting axis titles in 2D")
    OpenDatabase("../data/silo_%s_test_data/noise2d.silo"%SILO_MODE)
    AddPlot("Pseudocolor", "hardyglobal")
    DrawPlots()

    v = GetView2D()
    v.viewportCoords = (0.35, 0.95, 0.15, 0.95)
    SetView2D(v)

    a = GetAnnotationAttributes()
    a.axes2D.visible = 1
    a.axes2D.xAxis.label.visible = 0
    a.axes2D.yAxis.label.visible = 0
    a.axes2D.xAxis.title.visible = 1
    a.axes2D.yAxis.title.visible = 1
    a.xTitleFontHeight2D = 0.04
    a.yTitleFontHeight2D = 0.04
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_0_00")

    a.xAxisUserTitle2D = "New X Title"
    a.xAxisUserTitleFlag2D = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_0_01")

    a.yAxisUserTitle2D = "New Y Title"
    a.yAxisUserTitleFlag2D = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_0_02")

    a.xAxisUserUnits2D = "New X Units"
    a.xAxisUserUnitsFlag2D = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_0_03")

    a.yAxisUserUnits2D = "New Y Units"
    a.yAxisUserUnitsFlag2D = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_0_04")
    DeleteAllPlots()


#
# Test replacing 3D titles and units.
#
def Test3D():
    TestSection("Setting axis titles in 3D")
    OpenDatabase("../data/silo_%s_test_data/noise.silo"%SILO_MODE)
    AddPlot("Pseudocolor", "hardyglobal")
    DrawPlots()

    v = GetView3D()
    v.viewNormal = (-0.749133, -0.494511, 0.440747)
    v.focus = (0, 0, 0)
    v.viewUp = (-0.588718, 0.802033, -0.10077)
    v.viewAngle = 30
    v.parallelScale = 17.3205
    v.nearPlane = -34.641
    v.farPlane = 34.641
    v.imagePan = (0, 0)
    v.imageZoom = 1
    v.perspective = 1
    v.eyeAngle = 2
    v.centerOfRotationSet = 0
    v.centerOfRotation = (0, 0, 0)
    SetView3D(v)

    a = GetAnnotationAttributes()
    a.axes3D.visible = 1
    a.axes3D.xAxis.title.visible = 1
    a.axes3D.yAxis.title.visible = 1
    a.axes3D.zAxis.title.visible = 1
    a.xTitleFontHeight = 0.04 # 2x scale
    a.yTitleFontHeight = 0.04
    a.zTitleFontHeight = 0.04
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_00")

    a.xAxisUserTitle = "New X Title"
    a.xAxisUserTitleFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_01")

    a.yAxisUserTitle = "New Y Title"
    a.yAxisUserTitleFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_02")

    a.zAxisUserTitle = "New Z Title"
    a.zAxisUserTitleFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_03")

    a.xAxisUserUnits = "New X Units"
    a.xAxisUserUnitsFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_04")

    a.yAxisUserUnits = "New Y Units"
    a.yAxisUserUnitsFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_05")

    a.zAxisUserUnits = "New Z Units"
    a.zAxisUserUnitsFlag = 1
    SetAnnotationAttributes(a)
    SaveTestImage("axistitles_1_06")
    DeleteAllPlots()

def main():
    Test2D()
    Test3D()

# Run all of the tests
main()

# Exit the test
Exit()
