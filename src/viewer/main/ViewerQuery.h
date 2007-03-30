// ************************************************************************* //
//                              ViewerQuery.h                                //
// ************************************************************************* //

#ifndef VIEWER_QUERY_H
#define VIEWER_QUERY_H

#include <viewer_exports.h>
#include <SimpleObserver.h>
#include <ref_ptr.h>

// Forward declarations.
class AttributeSubject;
class Line;
class PlaneAttributes;
class PlotQueryInfo;
class SILRestrictionAttributes;
class ViewerPlot;
class ViewerWindow;
class avtToolInterface;
class avtVector;


// ****************************************************************************
//  Class: ViewerQuery
//
//  Purpose:  
//
//  Programmer: Kathleen Bonnell 
//  Creation:   June 10, 2002 
//
//  Modifications:
//    Kathleen Bonnell, Fri Jul 12 18:42:11 PDT 2002
//    Added width & height, for scaling purposes. Allow the
//    results window to be retrieved.
//
//    Kathleen Bonnell, Sat Jul 13 18:03:18 PDT 2002 
//    Added methods for handling tools. 
//
//    Kathleen Bonnell, Fri Jul 26 15:45:13 PDT 2002
//    Remove unused member origPlotQueryInfo.
//
//    Kathleen Bonnell, Mon Jul 29 09:36:35 PDT 2002  
//    Remove unnecessary methods InteractiveOn, InteractiveOff. 
//
//    Kathleen Bonnell, Thu Mar  6 15:15:30 PST 2003 
//    Added methods GetOriginatingWindow, GetOriginatingPlot, SendVisualCue, 
//    ReCreateLineout, UpdateLineFromSlice, Start/StopObservingPlot. 
//
// ****************************************************************************


class VIEWER_API ViewerQuery : public SimpleObserver
{
  public:
                     ViewerQuery(ViewerWindow *, ViewerWindow *, Line *);
                    ~ViewerQuery();

    bool             MatchResultsPlot(ViewerPlot *vp) const; 
    bool             MatchOriginatingPlot(ViewerPlot *vp) const; 

    bool             MatchResultsWindow(ViewerWindow *vw) const; 
    bool             MatchOriginatingWindow(ViewerWindow *vw) const; 

    void             DeleteOriginatingWindow();
    void             DeleteOriginatingPlot();
    void             DeleteVisualCue();

    virtual void     Update(Subject *) ;

    double           GetWidth() const;
    double           GetHeight() const;

    ViewerWindow    *GetResultsWindow() const;
    ViewerWindow    *GetOriginatingWindow() const;
    ViewerPlot      *GetOriginatingPlot() const;

    bool             CanHandleTool();
    bool             IsHandlingTool();
    bool             InitializeTool(avtToolInterface &ti);
    bool             HandleTool(const avtToolInterface &ti);
    void             DisableTool();
    void             SendVisualCue();
    void             ReCreateLineout();

    bool             UpdateLineFromSlice(PlaneAttributes *);

  private:
    void             CreateLineout();
    void             StartObservingPlot();
    void             StopObservingPlot();

    PlotQueryInfo   *resPlotQueryInfo;

    Line            *lineAtts;
    ViewerWindow    *originatingWindow;
    ViewerWindow    *resultsWindow;

    ViewerPlot      *resultsPlot;
    ViewerPlot      *originatingPlot;

    double           width;
    double           height;
    bool             handlingTool;

    PlaneAttributes *planeAtts;
};

typedef ref_ptr<ViewerQuery> ViewerQuery_p; 

#endif

