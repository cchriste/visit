# [VisIt Python Script]
#
# Driver to test example process id python expression.
#
# Usage:
#  visit -o path/to/rect2d.silo -nowin -cli -s visit_py_expr_procid_test.py
#

def script_path():
    # path magic to make sure we can locate the script file
    # no matter where we run visit from
    script_dir = os.path.split(__visit_source_file__)[0]
    return os.path.join(script_dir,"py_expr_procid.py")

def setup_save_win():
    swa= SaveWindowAttributes()
    swa.outputToCurrentDirectory = 1
    swa.outputDirectory = "."
    swa.fileName = "output_py_expr_procid_test"
    swa.family = 0
    swa.format = swa.PNG
    swa.width = 500
    swa.height = 500
    swa.screenCapture = 0
    swa.saveTiled = 0
    swa.quality = 100
    swa.progressive = 0
    swa.binary = 0
    swa.stereo = 0
    swa.compression = swa.PackBits
    swa.forceMerge = 0
    swa.resConstraint = swa.NoConstraint
    SetSaveWindowAttributes(swa)

if __name__ == "__main__":
    setup_save_win()
    DefinePythonExpression("myvar",["d"],file=script_path())
    AddPlot("Pseudocolor","myvar")
    DrawPlots()
    SaveWindow()
    sys.exit(0)

