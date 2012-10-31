# ----------------------------------------------------------------------------
#       Code to help in the VisIt test suite
#
#  Programmer: Jeremy Meredith
#  Date:       April 17, 2002
#
#  Modifications:
#    Hank Childs, Fri May 24 08:43:58 PDT 2002
#    Renamed SaveImageAtts to SaveWindowAtts.
#
#    Jeremy Meredith, Fri Jul 19 17:27:10 PDT 2002
#    Added python coloring code.  Added a third exit code for small errors.
#    Added code to write each individual test case to its own html file
#    and each test script to its own as well, and changed the formatting.
#
#    Jeremy Meredith, Thu Aug 29 15:10:45 PDT 2002
#    Improved the log file writing.
#
#    Jeremy Meredith, Fri Sep 13 17:11:48 PDT 2002
#    Made it brighten the difference image before saving to make it
#    easier to spot differences.  Made the per test case statistics
#    be in float format instead of decimal.
#
#    Hank Childs, Wed Nov 20 14:58:17 PST 2002
#    Explicitly test for divide by zero error.
#
#    Kathleen Bonnell, Fri Jun  6 12:09:41 PDT 2003  
#    Added TestText.
#
#    Jeremy Meredith, Mon Jun  9 17:53:36 PDT 2003
#    Moved the colorize-python code into its own module.
#    Added more advanced differencing to the textual comparisons.
#
#    Jeremy Meredith, Thu Jul 24 09:52:21 PDT 2003
#    Stopped saving baseline images.  It was messing things up.  Instead,
#    I made a default baseline image that clearly states "No Baseline Image".
#
#    Jeremy Meredith, Mon Aug 11 17:46:22 PDT 2003
#    Upped the quality level on output images to 90% (from the default of 75%)
#
#    Jeremy Meredith, Mon Aug 18 15:19:01 PDT 2003
#    Added timings.
#
#    Brad Whitlock, Mon Mar 22 13:51:17 PST 2004
#    Added TestSection.
#
#    Mark C. Miller, Tue Mar 30 15:48:25 PST 2004
#    Added new global, iactive, for interactive mode
#
#    Brad Whitlock, Tue Mar 30 15:56:38 PST 2004
#    I added code to create a graph of memory usage as the test runs.
#
#    Brad Whitlock, Fri Apr 2 10:00:09 PDT 2004
#    I fixed the memory tracking code so it should work when run in the
#    nightly test suite. Previously, it was failing because os.getlogin()
#    was throwing an exception, which prevented us from processing
#    any of the output from top. That resulted in empty memory plots.
#
#    Mark C. Miller, Tue May 11 20:21:24 PDT 2004
#    Changed scalable rendering controls to use activation mode
#
#    Jeremy Meredith, Thu Oct 21 13:24:51 PDT 2004
#    Difference images now are monochrome instead of grayscale.
#
#    Mark C. Miller, Mon Nov 29 18:52:41 PST 2004
#    Added code to do differences based upon check sums, if available
#    Made it so thumb and full size images for baseline and diff are NOT
#    generated in cases where test passes
#
#    Mark C. Miller, Mon Feb 14 20:24:23 PST 2005
#    Added code to deal with where to find ImageMagick convert utility
#    Removed code to re-define some CLI functions for HDF5 test mode
#
#    Mark C. Miller, Tue May 10 19:53:08 PDT 2005
#    Made it smarter about measuring differences
#
#    Mark C. Miller, Mon Jan 23 16:11:59 PST 2006
#    Changed default SaveWindowAttributes to NOT use screen capture. So,
#    had to explicitly invoke it here.
#
#    Mark C. Miller, Sat Feb 11 12:07:27 PST 2006
#    Force save in screen capture for background image
#
#    Mark C. Miller, Tue Nov 21 09:06:25 PST 2006
#    Changed code to remove userInfo from annotations to
#    SetDefaultAnnotationAttributes so it will take effect in all windows.
#    Re-organized code in this file to place functions at top and main 
#    execution lines at bottom
#
#    Mark C. Miller, Mon Nov 27 12:44:37 PST 2006
#    Work around bug in calling SetDefaultAnnotationAttributes by also
#    calling SetAnnotationAttributes
#
#    Mark C. Miller, hu Feb  8 17:13:14 PST 2007
#    Added population of 'silo' mode and logic to FilterTestText to
#    deal with test view used for silo's tests
#
#    Mark C. Miller, Tue Jan 29 16:37:53 PST 2008
#    Removed 'optimized' mode. Added -numdifftol command line option
#    and numdifftol global tolerance for relative numerical differences.
#
#    Tom Fogal, Wed Jan  6 18:06:09 MST 2010
#    Print out the import error so we have a chance of debugging it.
#
#    Jeremy Meredith, Tue Jan 19 11:43:41 EST 2010
#    Set the preferred plugin list to be exclusively Silo.  This mimics
#    the old behavior.
#
#    Mark C. Miller, Fri Jan 22 20:17:13 PST 2010
#    Added this comment to explain a previous update in which I added
#    externDbPaths variable and FindAndOpenDatabase() function.
#
#    Jeremy Meredith, Fri Mar 26 10:33:44 EDT 2010
#    It was decided that Silo should be a global preferred file format
#    for all users everywhere, so I removed the setting in this file.
#
#   Eric Brugger, Thu Apr 22 12:56:41 PDT 2010
#   I made several changes to the return code behavior of the script.  It
#   returns error code 119 if the test succeeded and the test had some skips.
#   It returns error code 120 if the test had acceptable differences and
#   had some skips.  It returns error code 113 if the differences were
#   unacceptable regardless of whether some tests were skipped.
#
# ----------------------------------------------------------------------------

import string
import sys
import time
import os
import glob
import subprocess
import thread

import HtmlDiff
import HtmlPython

from stat import *

pil_available = True
try:
    from PIL import Image, ImageChops, ImageStat
except ImportError, err:
    pil_available=False

# try to use visit_utils for logging?

###############################################################################
#
# Path helper Methods
#
###############################################################################

def abs_path(*args):
    """
    Helper for constructing absolute paths from a string, or lists of strings.
    """
    rargs = []
    for arg in args:
        if os.path.isabs(arg) or arg.count("/")  == 0:
            rargs.append(arg)
        else:
            toks = arg.split("/")
            rargs.extend(toks)
    res = pjoin(*rargs)
    return res

def out_path(*args):
    """
    Generates proper absolute path relative to the test suite results directory.
    """
    rargs = [TestEnv.visitResultDir]
    rargs.extend(args)
    return abs_path(*rargs)

def data_path(*args):
    """
    Generates proper absolute path relative to the 'data' directory.
    """
    rargs = [TestEnv.visitDataDir]
    rargs.extend(args)
    return abs_path(*rargs)

def silo_data_path(*args):
    """
    Helper that generates proper silo data absolute file paths.
    Incorporates SILO_MODE logic.
    """
    rargs = ["silo_%s_test_data" % SILO_MODE]
    rargs.extend(args)
    return data_path(*rargs)

def tests_path(*args):
    """
    Generates proper absolute path relative to the 'test/tests' directory.
    """
    rargs = [TestEnv.visitTopDir,"test","tests"]
    rargs.extend(args)
    return abs_path(*rargs)

def TestScriptPath():
    """
    Helper that provides the full path to the current script file

    Programmer: Cyrus Harrison
    Date:      Fri May 21 10:05:58 PDT 2010
    """
    script_file     = os.environ['VISIT_TEST_NAME']
    script_category = os.environ['VISIT_TEST_CATEGORY']
    script_dir      = tests_path(script_category,script_file)
    return script_dir

def sexe(cmd,ret_output=False,echo = False):
    """
    Helper for executing shell commands.
    """
    if echo:
        Log("[exe: %s]" % cmd)
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res =p.communicate()[0]
        return p.returncode,res
    else:
        return subprocess.call(cmd,shell=True)

def Log(msg):
    """
    Prints message to screen and also records to a log file.
    """
    print msg
    if (os.path.isfile("log")):
        log = open("log", 'a')
        log.write(msg)
        if msg.count("\n") == 0:
            log.write("\n")
        log.close()

# ----------------------------------------------------------------------------
# Function: GenFileNames
#
# Purpose:
#   Return the file names of the baseline, the current and the difference
#   files.
#
# Returns:
#   cur      The name of the current file.
#   diff     The name of the difference file.
#   base     The name of the baseline file.
#   altbase  The name of the mode specific baseline.
#   modeSpecific A flag indicating if the baseline is mode specific.
#
# Modifications:
#   Eric Brugger, Tue Apr 27 13:23:16 PDT 2010
#   I added the modeSpecific return value, which indicates if the baseline
#   image is mode specific.
#
# ----------------------------------------------------------------------------
def GenFileNames(file, ext):
    pyfilebase = TestEnv.pyfilebase 
    category   = TestEnv.category
    modeStr    = TestEnv.modeStr

    dcur_cat   = pjoin(TestEnv.visitResultDir,"current",category)
    dcur_base  = pjoin(dcur_cat,TestEnv.pyfilebase)
    ddiff_cat  = pjoin(TestEnv.visitResultDir,"diff",category)
    ddiff_base = pjoin(ddiff_cat,TestEnv.pyfilebase)

    for rdir in [dcur_cat, dcur_base, ddiff_cat, ddiff_base]:
        if not os.path.isdir(rdir):
            os.mkdir(rdir)
    # create file names
    cur  = pjoin(dcur_base,file + ext)
    diff = pjoin(ddiff_base,file + ext)
    base = pjoin(TestEnv.visitTestDir,"baseline",
                 category,pyfilebase,file + ext)
    altbase = ""
    modeSpecific = 0
    if modeStr != "":
        altbase = pjoin(TestEnv.visitTestDir,
                        "baseline",category,pyfilebase,modeStr,file + ext)
        if (os.path.isfile(altbase)):
            base = altbase
            modeSpecific = 1

    return (cur, diff, base, altbase, modeSpecific)

# ----------------------------------------------------------------------------
# Function: Test
#
# Purpose:
#   Write out the file, compare it to the baseline, thumbnail it,
#   and add it's data to the html
#
# Modifications:
#   Mark C. Miller, Mon Mar 29 19:37:26 PST 2004
#   Added alternate SaveWindowAttributes
#
#   Mark C. Miller, Tue Mar 30 15:48:25 PST 2004
#   Added pause for interacitve mode
#
#   Brad Whitlock, Tue Mar 30 16:38:52 PST 2004
#   Added code to sample memory.
#
#   Mark C. Miller, Mon Apr 12 16:34:50 PDT 2004
#   Added code to test against an alternate baseline if one exists
#
#   Jeremy Meredith, Tue May  4 13:26:41 PDT 2004
#   Catch exceptions from failing to open a baseline.  This can happen if
#   you make a clearcase element before putting an image into it.
#
#   Mark C. Miller, Tue Nov 28 23:13:08 PST 2006
#   Replaced maxerr, maxrms, maxpix with diffState. Added diff measure
#   indicating amount of diffs in pixels and not just count of diff pixels
#
#   Mark C. Miller, Tue Nov 28 23:50:15 PST 2006
#   Changed maxdiff to meddiff
#
#   Mark C. Miller, Wed Nov 29 08:19:52 PST 2006 
#   Changed meddiff to avgdiff
#
#   Sean Ahern, Thu Dec 20 14:48:14 EST 2007
#   Made diffState be a string so its easier to understand.
#
#   Eric Brugger, Thu Apr 22 12:56:41 PDT 2010
#   I made several changes to the return code behavior of the script.  It
#   returns error code 119 if the test succeeded and the test had some skips.
#   It returns error code 120 if the test had acceptable differences and
#   had some skips.  It returns error code 113 if the differences were
#   unacceptable regardless of whether some tests were skipped.
#
#   Eric Brugger, Tue Apr 27 13:23:16 PDT 2010
#   I enhanced the routine so that the text next to the large baseline image
#   indicates if it is a mode specific image or not.
#
#   Eric Brugger, Tue Jul 13 15:03:54 PDT 2010
#   I added the optional argument alreadySaved that indicates if the image
#   has already been saved.
#
# ----------------------------------------------------------------------------

def CalcDiffState(p_pixs, d_pixs, davg):
    if p_pixs != 0:
        dpix = d_pixs * 100.0 / p_pixs
        if dpix > TestEnv.pixdifftol:
            if davg > TestEnv.avgdifftol:
                diff_state = 'Unacceptable'
            else:
                diff_state = 'Acceptable'
        else:
            diff_state  = 'Acceptable'
    else:
        if d_pixs != 0:
            dpix = 1000000.0
            diff_state = 'Unacceptable'
        else:
            dpix = 0.0
            diff_state = 'None'
    return diff_state, dpix

def LogTestStart():
    # add test file info to log file
    msg  = "\n"
    msg += " - - - - - - - - - - - - - - -\n"
    msg += "  START:  Test script %s\n" % TestEnv.pyfilename
    msg += "\n"
    Log(msg)
    HTMLTestStart()

def LogTestExit(excode):
    msg  = "\n"
    msg += "  EXIT:   Test script %s\n" % TestEnv.pyfilename
    msg += "  EXCODE: %d\n" % excode
    msg += " - - - - - - - - - - - - - - -\n"
    msg += "\n"
    Log(msg)
    HTMLTestExit()

def HTMLTestStart():
    # set up our html output
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'wt')
    html.write("<SCRIPT TYPE=\"text/javascript\">\n")
    html.write("<!--\n")
    html.write("function popup(mylink, name)\n")
    html.write("{\n")
    html.write("if (! window.focus)return true;\n")
    html.write("var href;\n")
    html.write("if (typeof(mylink) == 'string')\n")
    html.write("   href=mylink;\n")
    html.write("else\n")
    html.write("   href=mylink.href;\n")
    html.write("window.open(href,name,'width=500,height=500,scrollbars=no');\n")
    html.write("return false;\n")
    html.write("}\n")
    html.write("//-->\n")
    html.write("</SCRIPT>\n")
    html.write("<html><head><title>Results for %s/%s</title></head>\n" % (TestEnv.category,TestEnv.pyfilename))
    html.write("<body bgcolor=\"#a0a0f0\">\n")
    html.write("<H1>Results of VisIt Regression Test - <a href=%s_%s_py.html>%s/%s</a></H1>\n" % (TestEnv.category,TestEnv.pyfilebase,TestEnv.category,TestEnv.pyfilename))
    html.write("<H2><a href=%s_%s_timings.html>(Full Timings)</a></H2>\n" % (TestEnv.category,TestEnv.pyfilebase))
    html.write("<table border>\n")
    html.write(" <tr>\n")
    html.write("  <td rowspan=2><b><i>Test Case</b></i></td>\n")
    html.write("  <td colspan=2 align=center><b><i>Errors</b></i></td>\n")
    html.write("  <td colspan=3 align=center><b><i>Images</b></i></td>\n")
    html.write(" </tr>\n")
    html.write(" <tr>\n")
    html.write("  <td>%Diffs</td>\n")
    html.write("  <td>Maximum</td>\n")
    html.write("  <td>Baseline</td>\n")
    html.write("  <td>Current</td>\n")
    html.write("  <td>Diff Map</td>\n")
    html.write(" </tr>\n")
    html.write("\n")

def HTMLTestExit():
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'a')
    html.write("</table>\n")
    html.write("</body>\n")
    html.write("</html>\n")


def LogTextTestResult(file,nchanges,nlines,failed,skip):
    # write data to the log file if there is one
    if failed:
        if skip:
            Log("    Test case '%s' SKIPPED" % file)
        else:
            Log("    Test case '%s' FAILED" % file)
    else:
        if nchanges < 0:
            Log("    Test case '%s' UNKNOWN" % file)
        else:
            Log("    Test case '%s' PASSED" % file)
    # write html result
    HTMLTextTestResult(file,nchanges,nlines,failed,skip)

def HTMLTextTestResult(file,nchanges,nlines,failed,skip):
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'a')
    # write to the html file
    color = "#00ff00"
    if failed:
        if skip:
            color = "#0000ff"
        else:
            color = "#ff0000"
    else:
        if (nchanges < 0):
            color = "#00ffff"
    html.write(" <tr>\n")
    html.write("  <td bgcolor=\"%s\"><a href=\"%s.html\">%s</a></td>\n" % (color, file, file))
    html.write("  <td colspan=5 align=center>%d modifications totalling %d lines</td>\n" % (nchanges,nlines))
    html.write(" </tr>\n")


def LogImageTestResult(file,diffState,modeSpecific,tPixs, pPixs, dPixs, dpix, davg):
    # write data to the log file if there is one
    if diffState == 'None':
        Log("    Test case '%s' PASSED" % file)
    elif diffState == 'Acceptable':
        Log("    Test case '%s' PASSED: #pix=%06d, #nonbg=%06d, #diff=%06d, ~%%diffs=%.3f, avgdiff=%3.3f" %
            (file, tPixs, pPixs, dPixs, dpix, davg))
    elif diffState == 'Unacceptable':
        Log("    Test case '%s' FAILED: #pix=%06d, #nonbg=%06d, #diff=%06d, ~%%diffs=%.3f, avgdiff=%3.3f" %
            (file, tPixs, pPixs, dPixs, dpix, davg))
    elif diffState == 'Skipped':
        Log("    Test case '%s' SKIPPED" % file)
    else:
        Log("    Test case '%s' UNKNOWN:#pix=UNK , #nonbg=UNK , #diff=UNK , ~%%diffs=UNK,  avgdiff=UNK")
    # write html result
    HTMLImageTestResult(file,diffState,modeSpecific,tPixs, pPixs,
                        dPixs, dpix, davg)

def Test(file, altSWA=0, alreadySaved=0):
    CheckInteractive()
    # for read only globals, we don't need to use "global"
    # we may need to use global for these guys
    #global maxds
    #global numskip

    (cur, diff, base, altbase, modeSpecific) = GenFileNames(file, ".png")

    # save the window in visit
    if alreadySaved == 0:
        if altSWA != 0:
            sa=altSWA
        else:
            sa=SaveWindowAttributes()
            sa.screenCapture=1
        sa.family   = 0
        sa.fileName = cur
        sa.format   = sa.PNG
        sa
        SetSaveWindowAttributes(sa)
        SaveWindow()

    diffState = 'Unknown'
    skip      = file in TestEnv.skipCases
    tPixs     = 0
    pPixs     = 0
    dPixs     = 0
    dpix      = 0.0
    davg      = 0.0
    if TestEnv.usePIL:
        (tPixs, pPixs, dPixs, davg) = DiffUsingPIL(file, cur, diff,
                                                   base, altbase)
        diffState, dpix = CalcDiffState(pPixs, dPixs, davg)
    if skip:
        diffState = 'Skipped'
        TestEnv.numskip += 1

    LogImageTestResult(file,diffState, modeSpecific,
                       dpix, tPixs, pPixs, dPixs, davg)

    # update maxmimum diff state
    diffVals = {
        'None' :         0,
        'Acceptable' :   1,
        'Unacceptable' : 2,
        'Unknown' :      3,
        'Skipped' :      0
    }
    TestEnv.maxds = max(TestEnv.maxds, diffVals[diffState])

# ----------------------------------------------------------------------------
# Function: HTMLImageTestResult
#
# Purpose:
#   Writes HTML stuff for a single test image 
#
# Modifications:
#   Mark C. Miller, Mon Jul  6 22:07:07 PDT 2009
#   Generate 'mouse-over' hrefs ONLY for case where there are diffs.
#   When there are no diffs, reference baseline instead of current. The
#   rationale for this later change is so that we can then create symlinks
#   for the html content instead of making copies and taking more disk space.
#
#   Eric Brugger, Tue Apr 27 13:23:16 PDT 2010
#   I added the modeSpecific argument, which causes the text next to the
#   baseline image to indicate if it is a mode specific image or not.
#
# ----------------------------------------------------------------------------

def HTMLImageTestResult(file,
                        diffState, modeSpecific,
                        dpix, tPixs, pPixs, dPixs, davg):
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'a')

    # write to the html file
    color = "#ffffff"
    if   diffState == 'None':           color = "#00ff00"
    elif diffState == 'Acceptable':     color = "#ffff00"
    elif diffState == 'Unacceptable':   color = "#ff0000"
    elif diffState == 'Skipped':        color = "#0000ff"
    else:                               color = "#ff00ff"
    html.write(" <tr>\n")
    html.write("  <td bgcolor=\"%s\"><a href=\"%s.html\">%s</a></td>\n" % (color, file, file))
    html.write("  <td align=center>%.2f</td>\n" % (dpix))
    html.write("  <td align=center>%.2f</td>\n" % (davg))
    if (diffState == 'Unknown'):
        html.write("  <td align=center>Not Available</td>\n")
        html.write("  <td align=center><a href=\"b_%s.png\" onclick='return popup(\"b_%s.png\",\"image\");'><img src=\"b_%s_thumb.png\"></a></td>\n" % (file,file,file))
        html.write("  <td align=center>Not Available</td>\n")
    elif (diffState != 'None'):
        html.write("  <td align=center><a href=\"b_%s.png\" onclick='return popup(\"b_%s.png\",\"image\");'><img src=\"b_%s_thumb.png\"></a></td>\n" % (file,file,file))
        html.write("  <td align=center><a href=\"c_%s.png\" onclick='return popup(\"c_%s.png\",\"image\");'><img src=\"c_%s_thumb.png\"></a></td>\n" % (file,file,file))
        html.write("  <td align=center><a href=\"d_%s.png\" onclick='return popup(\"d_%s.png\",\"image\");'><img src=\"d_%s_thumb.png\"></a></td>\n" % (file,file,file))
    else:
        html.write("  <td colspan=3 align=center><a href=\"b_%s.png\" onclick='return popup(\"b_%s.png\",\"image\");'><img src=\"b_%s_thumb.png\"></a></td>\n" % (file,file,file))
    html.write(" </tr>\n")

    # write the individual testcase
    testcase = open(out_path("html","%s.html" % file), 'wt')
    testcase.write("<html><head><title>Results for test case %s</title></head>\n" % file)
    testcase.write("<h1>Results for test case <i>%s</i></h1>\n" % file)
    testcase.write("<body bgcolor=\"#a080f0\">\n")
    testcase.write("<table border=5><tr><td></td></tr>\n")
    testcase.write("  <tr>\n")
    testcase.write("    <td align=center rowspan=9 bgcolor=%s>\n"%color)
    if (diffState == 'None' or diffState == 'Acceptable'):
        testcase.write("        <b><h1><i>Passed</i></h1></b>\n");
    elif (diffState == 'Unacceptable'):
        testcase.write("        <b><h1><i>Failed</i></h1></b>\n");
    elif (diffState == 'Skipped'):
        testcase.write("        <b><h1><i>Skipped</i></h1></b>\n");
    else:
        testcase.write("        <b><h1><i>Unknown</i></h1></b>\n");
    testcase.write("    </td>\n")
    if modeSpecific:
        testcase.write("    <td align=center>Mode<br>Specific<br>Baseline:</td>\n")
    else:
        testcase.write("    <td align=center>Baseline:</td>\n")
    if (diffState == 'None'):
        testcase.write("""    <td><img name="b" border=0 src="b_%s.png"></img></td>\n"""%file)
    elif (diffState == 'Unknown'):
        testcase.write("    <td>Not Available</td>\n")
    elif (diffState == 'Skipped'):
        testcase.write("""    <td><img name="b" border=0 src="b_%s.png"></img></td>\n"""%file)
    else:
        testcase.write("""    <td><a href="" onMouseOver="document.b.src='c_%s.png'" onMouseOut="document.b.src='b_%s.png'"><img name="b" border=0 src="b_%s.png"></img></a></td>\n"""%(file,file,file))
    testcase.write("  </tr>\n")
    testcase.write("  <tr>\n")
    testcase.write("    <td align=center>Current:</td>\n")
    if (diffState == 'None'):
        testcase.write("    <td>Same As Baseline</td>\n")
    elif (diffState == 'Unknown'):
        testcase.write("    <td>Not Available</td>\n")
    elif (diffState == 'Skipped'):
        testcase.write("    <td>Skipped</td>\n")
    else:
        testcase.write("""    <td><a href="" onMouseOver="document.c.src='b_%s.png'" onMouseOut="document.c.src='c_%s.png'"><img name="c" border=0 src="c_%s.png"></img></a></td>\n"""%(file,file,file))
    testcase.write("  </tr>\n")
    testcase.write("  <tr>\n")
    testcase.write("    <td align=center rowspan=7>Diff Map:</td>\n")
    if (diffState == 'None'):
        testcase.write("    <td rowspan=7>No Differences</td>\n")
    elif (diffState == 'Unknown'):
        testcase.write("    <td rowspan=7>Not Available</td>\n")
    elif (diffState == 'Skipped'):
        testcase.write("    <td rowspan=7>Skipped</td>\n")
    else:
        testcase.write("""    <td><a href="" onMouseOver="document.d.src='b_%s.png'" onMouseOut="document.d.src='d_%s.png'"><img name="d" border=0 src="d_%s.png"></img></a></td>\n"""%(file,file,file))
    testcase.write("    <td align=center><i>Error Metric</i></td>\n")
    testcase.write("    <td align=center><i>Value</i></td>\n")
    testcase.write("  </tr>\n")
    testcase.write("      <tr><td>Total Pixels</td>  <td align=right>%06d</td></tr>\n"%tPixs)
    testcase.write("      <tr><td>Non-Background</td><td align=right>%06d</td></tr>\n"%pPixs)
    testcase.write("      <tr><td>Different</td>     <td align=right>%06d</td></tr>\n"%dPixs)
    testcase.write("      <tr><td>%% Diff. Pixels</td><td align=right>%f</td></tr>\n"%dpix)
    testcase.write("      <tr><td>Avg. Diff</td><td align=right>%f</td></tr>\n"%davg)
    testcase.write("      <tr></tr>\n")
    testcase.write("  </tr>\n")
    testcase.write("</table>\n")
    testcase.write("</html>\n")
    testcase.close()

# ----------------------------------------------------------------------------
# Function: GetBackgroundImage 
#
# Purpose:
#   Returns the image of just VisIt's background without any plots
#
# ----------------------------------------------------------------------------

def GetBackgroundImage(file):

    notHiddenList = []
    activeList = []

    plots = ListPlots(1)
    plotInfos = string.split(plots,"#")
    for entry in plotInfos:

        if entry == "":
            continue;

        plotParts = string.split(entry,"|")
        plotHeader = plotParts[0]
        plotInfo = plotParts[1]

        # get CLI's plot id for this plot
        plotId = string.split(plotHeader,"[")
        plotId = plotId[1]
        plotId = string.split(plotId,"]")
        plotId = int(plotId[0])

        # get this plot's active & hidden status
        words = string.split(plotInfo,";")
        hidden = -1
        active = -1
        for word in words:
            if word == "hidden=0":
                hidden = 0
            elif word == "hidden=1":
                hidden = 1
            elif word == "active=0":
                active = 0
            elif word == "active=1":
                active = 1

        if active == 1:
            activeList.append(plotId)

        # if this plot is NOT hidden, hide it
        if hidden == 0:
            SetActivePlots((plotId))
            HideActivePlots()
            notHiddenList.append(plotId)

    # ok, all non-hidden plots have been hidden. So,
    # now save the background image
    oldSWA = SaveWindowAttributes()
    tmpSWA = SaveWindowAttributes()
    tmpSWA.family   = 0
    tmpSWA.fileName = out_path("current",file + "_bg.png")
    tmpSWA.format   = tmpSWA.PNG
    tmpSWA.screenCapture = 1
    SetSaveWindowAttributes(tmpSWA)
    SaveWindow()
    bkimage = Image.open(tmpSWA.fileName)

    # ok, now restore everything to the way it was
    # before we got here
    SetSaveWindowAttributes(oldSWA)
    SetActivePlots(tuple(notHiddenList))
    HideActivePlots()
    SetActivePlots(tuple(activeList))

    return bkimage

# ----------------------------------------------------------------------------
# Function: DiffUsingPIL 
#
# Purpose:
#   Diffs test results using PIL, outputs HTML, makes jpeg images,
#
# Modifications:
#   Jeremy Meredith, Tue Jun  7 12:14:11 PDT 2005
#   Fixed error reporting for missing baseline images.
#
#   Mark C. Miller, Mon Jul  6 22:07:07 PDT 2009
#   Modified to instead of always generating thumb of current (new) image
#   and only of baseline (old) and diffs when there are diffs to do the
#   opposite. That is always generate thumb of baseline (old) image and
#   current and diffs when there are diffs.
#
#   Mark C. Miller, Tue Jul 20 19:27:09 PDT 2010
#   Left in (commented out) line for color image stats. Will use in a later
#   update to compute max channel difference.
# ----------------------------------------------------------------------------

def DiffUsingPIL(file, cur, diff, baseline, altbase):

    # open it using PIL Image
    newimg = Image.open(cur)
    size = newimg.size;

    # open the baseline image
    try:
        if (os.path.isfile(altbase)):
            oldimg = Image.open(altbase)
            if (size != oldimg.size):
                Log("Error: Baseline and current images are different sizes... resizing baseline to compensate")
                oldimg = oldimg.resize(size, Image.BICUBIC)
        elif (os.path.isfile(baseline)):
            oldimg = Image.open(baseline)
            if (size != oldimg.size):
                Log("Error: Baseline and current images are different sizes... resizing baseline to compensate")
                oldimg = oldimg.resize(size, Image.BICUBIC)
        else:
            Log("Warning: No baseline image: %s" % baseline)
            oldimg = Image.open('nobaseline.pnm')
            oldimg = oldimg.resize(size, Image.BICUBIC)
    except:
        oldimg = Image.open('nobaseline.pnm')
        Log("Warning: Defective baseline image: %s" % baseline)
        oldimg = oldimg.resize(size, Image.BICUBIC)


    # create the difference image
    diffimg = ImageChops.difference(oldimg, newimg)
    #dstatc = ImageStat.Stat(diffimg) # stats of color image
    diffimg = diffimg.convert("L", (0.3333333, 0.3333333, 0.3333333, 0))

    # get some statistics
    dstat   = ImageStat.Stat(diffimg)
    dmin    = dstat.extrema[0][0]
    dmax    = dstat.extrema[0][1]
    dmean   = dstat.mean[0]
    dmedian = dstat.median[0]
    drms    = dstat.rms[0]
    dstddev = dstat.stddev[0]

    plotpixels = 0
    diffpixels = 0
    size = newimg.size
    totpixels = size[0] * size[1]

    mdiffimg = diffimg.copy()
    if (dmax > 0 and dmax != dmin):

        # brighten the difference image before we save it
        map = []
        map.append(0)
        for i in range(1,256): map.append(255)
        mdiffimg = mdiffimg.point(map)

        annotAtts = GetAnnotationAttributes()

        if (annotAtts.backgroundMode != 0 or 
            annotAtts.backgroundColor != (255, 255, 255, 255)):

            # we have to be really smart if we don't have a constant color
            # background
            backimg = GetBackgroundImage(file)

            # now, scan over pixels in oldimg counting how many non-background
            # pixels there are and how many diffs there are
            for col in range(0,size[0]):
                for row in range(0, size[1]):
                    newpixel = newimg.getpixel((col,row))
                    oldpixel = oldimg.getpixel((col,row))
                    backpixel = backimg.getpixel((col,row))
                    diffpixel = mdiffimg.getpixel((col,row))
                    if oldpixel != backpixel:
                        plotpixels = plotpixels + 1
                    if diffpixel == 255:
                        diffpixels = diffpixels + 1

        else:

            mdstat   = ImageStat.Stat(mdiffimg)
            oldimgm = oldimg.convert("L", (0.3333333, 0.3333333, 0.3333333, 0))
            map1 = []
            for i in range(0,254): map1.append(255)
            map1.append(0)
            map1.append(0)
            oldimgm = oldimgm.point(map1)
            mbstat   = ImageStat.Stat(oldimgm)
            diffpixels = int(mdstat.sum[0]/255)
            plotpixels = int(mbstat.sum[0]/255)

    mdiffimg.save(diff)

    # thumbnail size (w,h)
    thumbsize = (100,100)

    # create thumbnails and save jpegs
    oldthumb = oldimg.resize(   thumbsize, Image.BICUBIC)
    oldthumb.save(out_path("html","b_%s_thumb.png"%file));
    oldimg.save(out_path("html","b_%s.png"%file));
    if (dmax != 0):
        newthumb    = newimg.resize(   thumbsize, Image.BICUBIC)
        diffthumb   = mdiffimg.resize(  thumbsize, Image.BICUBIC)
        newthumb.save(out_path("html","c_%s_thumb.png"%file));
        diffthumb.save(out_path("html","d_%s_thumb.png"%file));
        newimg.save(out_path("html","c_%s.png"%file));
        mdiffimg.save(out_path("html","d_%s.png"%file));

    return (totpixels, plotpixels, diffpixels, dmean)

# ----------------------------------------------------------------------------
# Function: FilterTestText
#
# Purpose:
#   Filters words from the test text before it gets saved.
#
# Modifications:
#   Mark C. Miller, Tue Jan 29 18:57:45 PST 2008
#   Moved code to filter VISIT_TOP_DIR to top of routine to ensure it is
#   executed always, not just in the case where numdifftol is zero. I also
#   fixed cases where a floating point number occuring at the end of a
#   sentence ending in a period ('.') was not getting correctly interpreted
#   and filtered.
#
#   Mark C. Miller, Tue Jan 29 19:37:54 PST 2008
#   Adjusted case with the absolute value of the base values are very near
#   zero (e.g. less than square of numdifftol), to switch to absolute
#   diffs.
#
#   Mark C. Miller, Thu Jan 31 17:51:21 PST 2008
#   Modified the algorithm in a subtle way. Since the string.replace() calls
#   are applied repeatedly over the entire input string, it was possible for
#   an earlier replace to be corrupted by a later replace. The new algorithm
#   uses a two-step replace process. As numbers are found in the input string,
#   they are compared for magnitude of numerial difference to their counter
#   parts in the baseline string. If the difference is BELOW threshold, we aim
#   to replace the 'current' numerical value in inText with the baseline value
#   in baseText. This has the effect of preventing the numerical difference
#   from causing a REAL difference when the two strings are diff'd later on.
#   If the difference is NOT below threshold, we skip this replacement. That
#   has the effect of causing a REAL difference when the two strings are
#   diff'd later. When the replacement is performed (e.g. the numerical 
#   difference is below threshold), we perform the replacement in two steps.
#   In the first pass over the string, we replace each current value with
#   a unique replacement 'tag.' The string we use must be unique over all
#   words in inText. In the second pass, we replace each of these replacement
#   tags with the actual baseline string thereby making that word in the
#   string identical to the baseline result and effectively eliminating it
#   from effecting the text difference.
#
#   Mark C. Miller, Tue Mar  4 18:35:45 PST 2008
#   Fixed some issues with the replace algorithm. Changed how near-zero
#   diffs are handled back to 'ordinary' relative diff. Made it more graceful
#   if it is unable to import PIL. Made text diff'ing proceed without PIL.
#
#   Mark C. Miller, Tue Mar  4 19:53:19 PST 2008
#   Discovered that string managment was taking a non-trivial portion of
#   total test time for text-oriented tests. So, found a better way to handle
#   the replacements using string slicing. Now, replacements are handled as
#   we march word-for-word through the strings.
#
#   Mark C. Miller, Thu Mar  6 09:39:43 PST 2008
#   Made logic for relative diff clearer. Switched to use min operation in
#   denominator and switch order of min/abs there.
#
#   Mark C. Miller, Tue Jun  9 09:23:30 PDT 2009
#   Removed refs to old Clearcase VOB paths. Fixed setting of 'inStart' for
#   string replacement to use inWords[w] for search rather than inWordT
#   which is potentially altered due to translate call and may not be found
#   in tmpText.
# ----------------------------------------------------------------------------

def FilterTestText(inText, baseText):
    #
    # We have to filter out the absolute path information we might see in
    # this string. runtest passes the value for visitTopDir here.
    #
    inText = string.replace(inText, TestEnv.visitDataDir, "VISIT_TOP_DIR/data")
    inText = string.replace(inText, TestEnv.visitTestDir, "VISIT_TOP_DIR/test")
    numdifftol = TestEnv.numdifftol
    #
    # Only consider doing any string substitution if numerical diff threshold
    # is non-zero
    #
    if numdifftol != 0.0:

        tmpText = inText

        #
        # Break the strings into words. Pass over words looking for words that
        # form numbers. Whenever we have numbers, compute their difference
        # and compare it to threshold. If its above threshold, do nothing.
        # The strings will wind up triggering a text difference. If its below
        # threshold, eliminate the word from effecting text difference by
        # setting it identical to corresponding baseline word.
        #
        baseWords = string.split(baseText)
        inWords = string.split(tmpText)
        outText=""
        transTab = string.maketrans(string.digits, string.digits)
        inStart = 0
        for w in range(len(baseWords)):
            try:
                inWordT = string.translate(inWords[w], transTab, '><,()')
                baseWordT = string.translate(baseWords[w], transTab, '><,()')
                if inWordT.count(".") == 2 and inWordT.endswith(".") or \
                   baseWordT.count(".") == 2 and baseWordT.endswith("."):
                    inWordT = inWordT.rstrip(".")
                    baseWordT = baseWordT.rstrip(".")
                inStart = string.find(tmpText, inWords[w], inStart)

                #
                # Attempt to convert this word to a number. Exception indicates
                # it wasn't a number and we can move on to next word
                #
                inVal = string.atof(inWordT)
                baseVal = string.atof(baseWordT)

                #
                # Compute a relative difference measure for these two numbers
                # This logic was taken from www.math.utah.edu/~beebe/software/ndiff
                #
                if inVal == baseVal:
                    valDiff = 0
                elif inVal == 0 and baseVal != 0:
                    valDiff = numdifftol # treat as above threshold 
                elif inVal != 0 and baseVal == 0:
                    valDiff = numdifftol # treat as above threshold 
                else:
                    valDiff = abs(inVal - baseVal) / min(abs(inVal), abs(baseVal))

                #
                # We want to ignore diffs that are deemed below threshold given
                # the relative diff. measure above. To affect this, we need to
                # replace the numbers in the input text that differ with their
                # cooresponding numbers in the baseline text. This will have the
                # effect of making the HTML difference ignore this value.
                # So, we do this replace only if the diff is non-zero and less
                # than threshold.
                #
                if valDiff > 0 and valDiff < numdifftol:
                    tmpText = tmpText[:inStart] + baseWordT + tmpText[inStart+len(inWordT):]

                inStart = inStart + len(inWordT)

            #
            # This word wasn't a number, move on
            #
            except ValueError:
                # ignore exceptions
                pass

        return tmpText

    else:

        return inText

def CheckInteractive():
    """
    Helper which pauses if we are in interactive mode.
    """
    # if interactive, pause for user
    if TestEnv.iactive:
        print "***********************"
        print "***********************"
        print "***********************"
        print "Saving %s"%file
        print "Hit Any Key To Continue"
        print "***********************"
        print "***********************"
        print "***********************"
        next = sys.stdin.read(1)

# ----------------------------------------------------------------------------
# Function: TestText
#
# Purpose:
#   Write out text to file, diff it with the baseline, and add it's data to 
#   the html
#
# Modifications:
#   Brad Whitlock, Tue Mar 30 16:39:43 PST 2004
#   Added code to sample memory.
#
#   Mark C. Miller, Tue May 25 14:29:40 PDT 2004
#   Added code to support interactive mode
#
#   Hank Childs, Sun Mar 27 14:34:30 PST 2005
#   Fix typo with testing mode specific baselines.
#
#   Jeremy Meredith, Tue Jun  7 12:09:12 PDT 2005
#   Added support for missing baseline text files again.
#
#   Brad Whitlock, Mon Nov 21 13:41:19 PST 2005
#   I made sure that it uses the mode-specific baseline if one exists.
#
#   Eric Brugger, Thu Apr 22 12:56:41 PDT 2010
#   I made several changes to the return code behavior of the script.  It
#   returns error code 119 if the test succeeded and the test had some skips.
#   It returns error code 120 if the test had acceptable differences and
#   had some skips.  It returns error code 113 if the differences were
#   unacceptable regardless of whether some tests were skipped.
#
#   Eric Brugger, Tue Apr 27 13:23:16 PDT 2010
#   I enhanced the routine so that the text next to the large baseline image
#   indicates if it is a mode specific image or not.
#
# ----------------------------------------------------------------------------
def TestText(file, inText):
    CheckInteractive()
    global maxds
    global numskip

    # create file names
    (cur, diff, base, altbase, modeSpecific) = GenFileNames(file, ".txt")

    if os.path.isfile(base):
        fin = open(base)
        baseText = fin.read()
    else:
        Log("Warning: No baseline text file: %s" % base)
        base = "notext.txt"
        baseText = "notext"

    # Filter out unwanted text
    inText = FilterTestText(inText, baseText)

    # save the current text output 
    fout = open(cur, 'w')
    fout.write(inText)
    fout.close()

    nchanges = 0
    nlines   = 0

    # diff the baseline and current text files
    d = HtmlDiff.Differencer(base, cur)
    # change to use difflib
    (nchanges, nlines) = d.Difference(out_path("html","%s.html"%file), file)

    # save the diff output 
    # TODO, this wont work on WINDOWS!
    # we can use difflib
    diff_cmd = "diff " + base + " " + cur
    r,diff_out = sexe(diff_cmd,ret_output = True)
    fout = open(diff, 'w')
    fout.write(diff_out)
    fout.close()

    # did the test fail?
    failed = (nchanges > 0)
    skip   = file in TestEnv.skipCases

    LogTextTestResult(file,nchanges,nlines,failed,skip)

    # Increment the number of skips if appropriate
    if skip:
        TestEnv.numskip += 1

    # set error codes
    # TODO: This doesn't seem correct, these should accumulate
    if failed:
        if skip:
            TestEnv.maxds = 0
        else:
            TestEnv.maxds = 2

# ----------------------------------------------------------------------------
# Function: TestSection
#
# Purpose:
#   Write a section header into the results table so it is easier to understand
#   the results for a large test.
#
# ----------------------------------------------------------------------------

def TestSection(sectionName):
    LogSectionStart(sectionName)

def LogSectionStart(sectionName):
    Log("    BEGIN SECTION: %s" % sectionName)
    HTMLSectionStart(sectionName)

def HTMLSectionStart(sectionName):
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'a')
    html.write(" <tr>\n")
    html.write("  <td colspan=6 align=center bgcolor=\"#0000ff\"><font color=\"#ffffff\"><b>%s</b></font></td>\n" % sectionName)
    html.write(" </tr>\n")


# ----------------------------------------------------------------------------
# Function: Exit
#
# Purpose:
#   Exit with the appropriate error code.  Must be called at end of test cases.
#
# Modifications:
#   Eric Brugger, Thu Apr 22 12:56:41 PDT 2010
#   I made several changes to the return code behavior of the script.  It
#   returns error code 119 if the test succeeded and the test had some skips.
#   It returns error code 120 if the test had acceptable differences and
#   had some skips.  It returns error code 113 if the differences were
#   unacceptable regardless of whether some tests were skipped.
# ----------------------------------------------------------------------------

def Exit(excode=0):
    html = open(out_path("html","%s_%s.html" % (TestEnv.category, TestEnv.pyfilebase)), 'a')
    # future mem tracking logic will need to clear engine cache:
    #ClearCacheForAllEngines()
    rcode = None
    if TestEnv.iactive == 0:
        if (excode):
            rcode = excode
    if rcode is None and TestEnv.maxds == 0:
        if TestEnv.numskip == 0:
            rcode = 111
        else:
            rcode = 119
    if rcode is None and TestEnv.maxds == 1:
        if TestEnv.numskip == 0:
            rcode = 112
        else:
            rcode = 120
    if rcode is None and TestEnv.maxds == 2:
        rcode = 113
    if rcode is None:
        rcode  = 114
    LogTestExit(rcode)
    sys.exit(rcode)

def TurnOnAllAnnotations(givenAtts=0):
    """
    Turns on all annotations.

    Either from the default instance of AnnotationAttributes,
    or using 'givenAtts'.
    """
    if (givenAtts == 0):
        a = AnnotationAttributes()
    else:
        a = givenAtts
    a.axes2D.visible = 1
    a.axes3D.visible = 1
    a.axes3D.triadFlag = 1
    a.axes3D.bboxFlag = 1
    a.userInfoFlag = 0
    a.databaseInfoFlag = 1
    a.legendInfoFlag = 1
    SetAnnotationAttributes(a)

def TurnOffAllAnnotations(givenAtts=0):
    """
    Turns off all annotations.

    Either from the default instance of AnnotationAttributes,
    or using 'givenAtts'.
    """
    if (givenAtts == 0):
        a = AnnotationAttributes()
    else:
        a = givenAtts
    a.axes2D.visible = 0
    a.axes3D.visible = 0
    a.axes3D.triadFlag = 0
    a.axes3D.bboxFlag = 0
    a.userInfoFlag = 0
    a.databaseInfoFlag = 0
    a.legendInfoFlag = 0
    SetAnnotationAttributes(a)


def FindAndOpenDatabase(dbname, extraPaths=()):
    """
    Searches in places where we're likely to find data files that we do NOT
    maintain in VisIt's repo. One needs to use the FindAndOpenDatabase()
    method in place of the OpenDatabase() method to invoke the behavior
    of searching in these dirs for the database file to open.
    """
    externalDbPaths=(data_path(),
                     "/project/projectdirs/visit/data",
                      "/usr/gapps/visit/data",
                      "/home/visit/data")
    for p in externalDbPaths + extraPaths:
        abs_dbname = "%s/%s"%(p,dbname)
        if os.path.isfile(abs_dbname):
            return OpenDatabase(abs_dbname), abs_dbname
    Log("Unable to OpenDatabase \"%s\" at any of the specified paths.\n" % dbname)
    return 0, ""

#############################################################################
#   Argument/Environment Processing
#############################################################################

class TestEnv(object):
    """
    Class that holds high level environment options.

    Replaces old use of old global vars.
    """
    # interactive?
    iactive     = 0
    # pil related
    usePIL      = 1
    avgdifftol  = 0.0
    pixdifftol  = 0
    numdifftol  = 0.0
    # mode inof 
    serial    = 1
    scalable  = 0
    parallel  = 0
    silo      = 0
    modeStr   = ""
    SILO_MODE = "hdf5"
    # these are outputs
    maxds       = 0
    numskip     = 0
    @classmethod
    def Setup(cls):
        # Process command line arguments.
        for arg in sys.argv:
            if (arg == "-noPIL"):
                cls.usePIL = 0
            if (arg == "-interactive"):
                cls.iactive = 1
            subargs = string.split(arg,"=")
            if (subargs[0] == "-pixdiff"):
                cls.pixdifftol = float(subargs[1])
            if (subargs[0] == "-avgdiff"):
                cls.avgdifftol = float(subargs[1])
            if (subargs[0] == "-numdiff"):
                cls.numdifftol = float(subargs[1])
        cls.visitTopDir    = os.environ['VISIT_TOP_DIR']
        cls.visitDataDir   = os.path.abspath(pjoin(cls.visitTopDir,"data"))
        cls.visitTestDir   = os.path.abspath(pjoin(cls.visitTopDir,"test"))
        cls.visitResultDir = cls.visitTestDir
        cls.pyfilename     = os.environ['VISIT_TEST_NAME']
        cls.pyfilebase     = cls.pyfilename[:-3]
        cls.category       = os.environ['VISIT_TEST_CATEGORY']
        cls.modes          = string.split(os.environ['VISIT_TEST_MODES'],",")
        cls.skipCases      = string.split(os.environ['VISIT_TEST_SKIP_CASES'],",")
        if os.environ.has_key("VISIT_RESULT_DIR"):
            cls.visitResultDir = os.path.abspath(os.environ["VISIT_RESULT_DIR"])
        # parse modes for various possible modes
        for mode in cls.modes:
            if cls.modeStr == "":
                cls.modeStr = mode
            else:
                if mode != "":
                    cls.modeStr = cls.modeStr + "_" + mode
            if mode == "scalable":
                cls.scalable = 1
            if mode == "parallel":
                cls.parallel = 1
                cls.serial = 0
            if mode == "silo":
                cls.silo = 1
            if mode == "pdb":
                cls.SILO_MODE = "pdb"
            if cls.usePIL and not pil_available:
                Log("WARNING: unable to import modules from PIL: %s" % str(err))
                cls.usePIL = 0

def InitTestEnv():
    """
    Sets up VisIt for a test.
    """
    # main setup
    TestEnv.Setup()
    # remove the user name
    annot = AnnotationAttributes()
    annot.userInfoFlag = 0
    SetDefaultAnnotationAttributes(annot)
    SetAnnotationAttributes(annot)
    # set scalable rendering mode if desired
    if TestEnv.scalable:
        ra = GetRenderingAttributes()
        ra.scalableActivationMode = ra.Always
        SetRenderingAttributes(ra)
    else:
        ra = GetRenderingAttributes()
        ra.scalableActivationMode = ra.Never
        SetRenderingAttributes(ra)
    # start parallel engine if parallel
    haveParallelEngine = 1L
    if TestEnv.parallel:
        haveParallelEngine = OpenComputeEngine("localhost", ("-np", "2"))
    if haveParallelEngine == 0L:
        Exit()
    else:
        OpenComputeEngine("localhost")
    # Automatically turn off all annotations
    # This is to prevent new tests getting committed that
    # are unnecessarily dependent on annotations.
    TurnOffAllAnnotations()
    # make sure the html dir exists (for compat w/ old runtest)
    if not os.path.isdir(out_path("html")):
        os.mkdir(out_path("html"))
    # colorize the source file, and write to an html file
    HtmlPython.ColorizePython(TestEnv.visitTestDir,
                              TestEnv.visitResultDir,
                              TestEnv.category,
                              TestEnv.pyfilename,
                              TestEnv.pyfilebase)
    LogTestStart()

# keep this as a global for now.
SILO_MODE = TestEnv.SILO_MODE

InitTestEnv()





