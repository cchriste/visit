# ----------------------------------------------------------------------------
#  CLASSES: nightly
#
#  Test Case:  saveformats.py
#
#  Tests:      save window formats 
#
#  Programmer: Mark C. Miller 
#  Date:       September 20, 2005 
#
#  Modifications:
#    Jeremy Meredith, Mon Apr 23 14:07:35 EDT 2007
#    Don't create a mesh plot when we're saving geometry formats, since
#    the new behavior is to separate plots into different files (as they
#    typically have different variables and geometry types).
#
#    Mark C. Miller, Wed Jan 20 07:37:11 PST 2010
#    Added ability to swtich between Silo's HDF5 and PDB data.
#
#    Cyrus Harrison, Tue Feb  2 10:55:43 PST 2010
#    Fixed problem w/ setting active window that allowed errors to propagate
#    between test cases.
#
# ----------------------------------------------------------------------------
swa=SaveWindowAttributes()
swa.family = 0

AddWindow()
SetActiveWindow(1)

# I=Image, G=Geometry, C=Curve formats
CFormats=[swa.CURVE, swa.POSTSCRIPT, swa.ULTRA]
IFormats=[swa.BMP, swa.JPEG, swa.PNG, swa.PPM, swa.RGB, swa.TIFF]
GFormats=[swa.STL, swa.OBJ, swa.VTK]

a = AnnotationAttributes()
a.userInfoFlag = 0
a.databaseInfoFlag = 0
a.legendInfoFlag = 0
SetAnnotationAttributes(a)

def TestSaveFormat(fmt):
    SetActiveWindow(1)
    mode = ""
    result = "Failed\n"
    (ext, isC, isI, isG) = SaveFileInfo(fmt)
    swatmp = swa
    swatmp.format = fmt
    if isI:
        swatmp.fileName = "current/saveformat_tmp.%s"%ext
    else:
        swatmp.fileName = "current/saveformat_tmp"
    SetSaveWindowAttributes(swatmp)
    SaveWindow()

    # depending on the type of format this is, try to
    # read the file we just created back into VisIt and
    # put up a plot in window 2. If that succeeds, we'll
    # say this format's save actually worked
    if isC:
        if FileExists(swatmp.fileName+"."+ext, 1, 0):
            if ext == "ps":
                result = "Passed\n" # can only test existence for ps
            else:
                SetActiveWindow(2)
                if OpenDatabase(swatmp.fileName+"."+ext):
                    AddPlot("Curve","curve")
                    if DrawPlots():
                        result = "Passed\n"
                    DeleteAllPlots()
                    CloseDatabase(swatmp.fileName+"."+ext)
    elif isI:
        if swatmp.screenCapture == 0:
            mode = "offscreen_"
        tiffFileName = "current/saveformat_tmp.tif"
        tiffFileExists = 0
        imageFileExists = FileExists(swatmp.fileName, 1, 0)
        if imageFileExists:
            os.system("%s %s -compress none %s"%(imgConverter, swatmp.fileName, tiffFileName))
            tiffFileExists = FileExists(tiffFileName, 1, 0)
        if tiffFileExists:
            SetActiveWindow(2)
            if OpenDatabase(tiffFileName):
                AddPlot("Pseudocolor","red")
                if DrawPlots():
                    result = "Passed\n"
                DeleteAllPlots()
                CloseDatabase(tiffFileName)
    elif isG:
        if FileExists(swatmp.fileName+"."+ext, 1, 0):
            if ext == "stl":
                meshName = "STL_mesh"
            elif ext == "obj":
                meshName = "OBJMesh"
            elif ext == "vtk":
                meshName = "mesh"
            SetActiveWindow(2)
            if OpenDatabase(swatmp.fileName+"."+ext):
                AddPlot("Mesh",meshName)
                if DrawPlots():
                    result = "Passed\n"
                DeleteAllPlots()
                CloseDatabase(swatmp.fileName+"."+ext)
    TestText("saveformat_%s%s"%(mode,ext), result)
    SetActiveWindow(1)

TestSection("Curve Formats")
OpenDatabase("../data/curve_test_data/c062.curve")
AddPlot("Curve", "going_down")
DrawPlots()
for f in CFormats:
    TestSaveFormat(f)

TestSection("Image Formats via Screen Capture")
DeleteAllPlots()
CloseDatabase("../data/curve_test_data/c062.curve")
OpenDatabase("../data/silo_%s_test_data/multi_rect2d.silo"%SILO_MODE)
AddPlot("Mesh", "mesh1")
AddPlot("Pseudocolor", "d")
DrawPlots()
slider = CreateAnnotationObject("TimeSlider")
for f in IFormats:
    TestSaveFormat(f)

TestSection("Image Formats via Off Screen")
swa.screenCapture = 0
for f in IFormats:
    TestSaveFormat(f)
swa.screenCapture = 1

TestSection("Geometry Formats")
DeleteAllPlots()
CloseDatabase("../data/silo_%s_test_data/multi_rect2d.silo"%SILO_MODE)
OpenDatabase("../data/silo_%s_test_data/globe.silo"%SILO_MODE)
AddPlot("Pseudocolor", "dx")
DrawPlots()
for f in GFormats:
    TestSaveFormat(f)

Exit()
