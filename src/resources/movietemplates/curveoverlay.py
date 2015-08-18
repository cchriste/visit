###############################################################################
# Function: Sequence1Frames_set_timeslider
#
# Purpose:
#   This is a callback function for sequence 1's IterateCallbackAndSaveFrames
#   function. This function sets the time and updates the time slider so
#   it has the right time value.
#
# Programmer: Brad Whitlock
# Creation:   Thu Nov 16 11:46:31 PDT 2006
#
# Modifications:
#   Brad Whitlock, Wed Apr 16 16:44:26 PDT 2008
#   Add the current time to the movie template so we eventually know all of
#   the times.
#
###############################################################################

def Sequence1Frames_set_timeslider(i, cbdata):
    ts = cbdata[0]
    ret = SetTimeSliderState(i)
    Query("Time")
    time = GetQueryOutputValue()
    ts.text = "Time = %1.5f" % time
    # Add the time to the movie template
    cmt = cbdata[1]
    cmt.databaseTimes = cmt.databaseTimes + [time]
    return ret

###############################################################################
# Function: Sequence1Frames_set_timeslider
#
# Purpose:
#   This is a callback function for sequence 2's IterateCallbackAndSaveFrames
#   function. This function lets us adjust the clip plane as a function of 
#   the number of time states and save out an image each time.
#
# Programmer: Brad Whitlock
# Creation:   Thu Nov 16 11:46:31 PDT 2006
#
# Modifications:
#   Brad Whitlock, Wed Apr 16 16:45:24 PDT 2008
#   Changed calculation of t so that it can be non-uniform.
#
#   Brad Whitlock, Thu Nov 11 15:21:21 PST 2010
#   We have to DrawPlots after setting operator options now.
#
###############################################################################

def Sequence2Frames_clip_cb(i, cbdata):
    nts = cbdata[0]
    clip = cbdata[1]
    xmin = cbdata[2]
    xmax = cbdata[3]
    vc = cbdata[4]
    allTimes = cbdata[5]
    dT = allTimes[-1] - allTimes[0]
    t = (allTimes[i] - allTimes[0]) / dT
    newX = t * (xmax - xmin) + xmin
    clip.plane1Origin = (newX, 0, 0)
    ret = SetOperatorOptions(clip)
    DrawPlots()
    SetViewCurve(vc)
    return ret

###############################################################################
# Class: OverlayCurveMovieTemplate
#
# Purpose:
#   This is movie template class creates a movie of a FilledBoundary plot
#   and a Curve plot that animates over time.
#
# Programmer: Brad Whitlock
# Creation:   Thu Nov 16 11:46:31 PDT 2006
#
# Modifications:
#   Brad Whitlock, Wed Apr 16 16:45:49 PDT 2008
#   Added databaseTimes member.
#
###############################################################################

class OverlayCurveMovieTemplate(VisItMovieTemplate):
    def __init__(self, mm, tr):
        super(OverlayCurveMovieTemplate, self).__init__(mm, tr)
        self.databaseTimes = []

    ###########################################################################
    # Function: Sequence1Frames
    #
    # Purpose:
    #   This method creates the frames for sequence 1.
    #
    # Programmer: Brad Whitlock
    # Creation:   Thu Nov 16 11:46:31 PDT 2006
    #
    # Modifications:
    #   Brad Whitlock, Wed Apr 16 16:46:27 PDT 2008
    #   Pass self to the sequence1 callback.
    #
    #   Brad Whitlock, Tue Apr 22 15:52:52 PDT 2008
    #   Turn off the query output.
    #
    ###########################################################################

    def Sequence1Frames(self, formats, percents):
        self.Debug(1, "OverlayCurveMovieTemplate.Sequence1Frames: begin")
        options = self.sequence_data["SEQUENCE_1"]

        # Set up the plots.
        DeleteAllPlots()
        OpenDatabase(options["DATABASE"])
        if options["PLOT_TYPE"] == 0:
            if AddPlot("FilledBoundary", options["PLOT_VAR"]) == 0:
                raise self.error("The FilledBoundary plot could not be created for "
                            "sequence 1.")
        else:
            if AddPlot("Pseudocolor", options["PLOT_VAR"]) == 0:
                raise self.error("The Pseudocolor plot could not be created for "
                            "sequence 1.")
        DrawPlots()

        # Set the background color.
        annot = GetAnnotationAttributes()
        annot.foregroundColor = (255, 255, 255, 255)
        annot.gradientColor1 = options["GRADIENT_BGCOLOR1"]
        annot.gradientColor2 = options["GRADIENT_BGCOLOR2"]
        annot.gradientBackgroundStyle = annot.TopToBottom
        annot.backgroundMode = annot.Gradient
        # Turn off certain annotations.
        annot.userInfoFlag = 0
        annot.databaseInfoFlag = 0
        annot.legendInfoFlag = 0
        # Set the axis names
        annot.axes2D.xAxis.title.title = options["XAXIS_TEXT"]
        annot.axes2D.yAxis.title.title = options["YAXIS_TEXT"]
        annot.axes2D.xAxis.title.userTitle = 1
        annot.axes2D.yAxis.title.userTitle = 1
        SetAnnotationAttributes(annot)

        # Change the viewport
        v = GetView2D()
        v.viewportCoords = (0.1, 0.95, 0.35, 0.95)
        SetView2D(v)

        ts = CreateAnnotationObject("TimeSlider")

        classification = CreateAnnotationObject("Text2D")
        classification.text = options["CLASSIFICATION_TEXT"]
        classification.useForegroundForTextColor = 0
        classification.textColor = options["CLASSIFICATION_TEXTCOLOR"]
        classification.position = (0.80, 0.97)
        classification.height = 0.02
        classification.fontBold = 1

        title = CreateAnnotationObject("Text2D")
        title.text = options["TITLE"]
        title.position = (0.01, 0.955)
        title.height = 0.03
        title.fontBold = 1

        # Save the frames.
        SuppressQueryOutputOn()
        cb_data = (TimeSliderGetNStates(), Sequence1Frames_set_timeslider, (ts, self))
        ret = self.IterateCallbackAndSaveFrames(cb_data, "seq1", formats, percents, "Generating sequence 1 frames")
        SuppressQueryOutputOff()

        DeleteAllPlots()
        ts.Delete()
        classification.Delete()
        title.Delete()

        self.Debug(1, "OverlayCurveMovieTemplate.Sequence1Frames: end")
        return (ret, "seq1", GetAnnotationAttributes().backgroundColor)

    ###########################################################################
    # Function: Sequence2Frames
    #
    # Purpose:
    #   This method creates the frames for sequence 2.
    #
    # Programmer: Brad Whitlock
    # Creation:   Thu Nov 16 11:46:31 PDT 2006
    #
    # Modifications:
    #   Brad Whitlock, Wed Apr 16 16:47:11 PDT 2008
    #   Pass the databaseTimes to the sequence 2 callback.
    #
    #   Brad Whitlock, Thu Nov 11 15:25:51 PST 2010
    #   Fix annotation function calls that were deprecated.
    #
    ###########################################################################

    def Sequence2Frames(self, formats, percents):
        self.Debug(1, "OverlayCurveMovieTemplate.Sequence2Frames: begin")
        options = self.sequence_data["SEQUENCE_2"]

        # Determine the number of time steps in the first sequence's database.
        options1 = self.sequence_data["SEQUENCE_1"]
        OpenDatabase(options1["DATABASE"])
        nts = TimeSliderGetNStates()
        CloseDatabase(options1["DATABASE"])
        DeleteAllPlots()
        self.DeleteAllAnnotationObjects()

        # Set up the Curve plot.
        OpenDatabase(options["CURVE_DATABASE"])
        AddPlot("Curve", options["CURVE_VARIABLE"])
        DrawPlots()
        cAtts = CurveAttributes(1)
        cAtts.showLabels = 0
        SetPlotOptions(cAtts)
        ResetView()
        vc = GetViewCurve()
        vc.viewportCoords = (0.1, 0.95, 0.15, 1.)

        # Get the Curve plot extents
        Query("SpatialExtents")
        extents = GetQueryOutputValue()
        self.Debug(5, "extents=" + str(extents))
        AddOperator("Clip")
        clip = ClipAttributes()
        clip.funcType = clip.Plane
        clip.plane1Status = 1
        clip.plane2Status = 0
        clip.plane3Status = 0
        clip.plane1Origin = (extents[0], 0, 0)
        clip.plane1Normal = (1, 0, 0)
        clip.planeInverse = 0
        SetOperatorOptions(clip)
        DrawPlots()

        # Set the background color.
        annot = GetAnnotationAttributes()
        annot.backgroundMode = annot.Solid
        annot.foregroundColor = (255, 255, 255, 255)
        annot.backgroundColor = (0, 0, 0, 255)
        # Turn off most annotations.
        annot.userInfoFlag = 0
        annot.databaseInfoFlag = 0
        annot.legendInfoFlag = 0
        annot.axes2D.xAxis.title.visible = 0
        annot.axes2D.yAxis.title.visible = 0
        annot.axes2D.xAxis.label.visible = 0
        annot.axes2D.yAxis.label.visible = 0
        SetAnnotationAttributes(annot)

        title = CreateAnnotationObject("Text2D")
        title.text = options["CURVE_TITLE"]
        title.position = (0.11, 0.88)
        title.height = 0.1
        title.fontBold = 1

        # Save the frames. This will be done by some other thing so the 
        # will have the viewport names worked in.
        cb_data = (nts, Sequence2Frames_clip_cb, (nts, clip, extents[0], extents[1], vc, self.databaseTimes))
        ret = self.IterateCallbackAndSaveFrames(cb_data, "seq2", formats, percents, "Generating sequence 2 frames")

        title.Delete()
        DeleteAllPlots()

        # Reset the database times
        self.databaseTimes = []

        self.Debug(1, "OverlayCurveMovieTemplate.Sequence2Frames: end")
        return (ret, "seq2", GetAnnotationAttributes().backgroundColor)

    ###########################################################################
    # Function: HandleScriptingSequence
    #
    # Purpose:
    #   This method invokes the appropriate routine for creating sequence
    #   frames.
    #
    # Programmer: Brad Whitlock
    # Creation:   Thu Nov 16 11:46:31 PDT 2006
    #
    # Modifications:
    #
    ###########################################################################

    def HandleScriptingSequence(self, seqName, formats, percents):
        ret = 0
        if seqName == "SEQUENCE_1":
            ret = self.Sequence1Frames(formats, percents)
        elif seqName == "SEQUENCE_2":
            ret = self.Sequence2Frames(formats, percents)
        return ret

# Public
def InstantiateMovieTemplate(moviemaker, templateReader):
    return OverlayCurveMovieTemplate(moviemaker, templateReader)
