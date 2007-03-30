// ************************************************************************* //
//                             ViewerPlotFactory.h                           //
// ************************************************************************* //

#ifndef VIEWER_PLOT_FACTORY_H
#define VIEWER_PLOT_FACTORY_H
#include <viewer_exports.h>
#include <avtSILRestriction.h>

class AttributeSubject;
class ViewerPlot;
class ViewerPlotPluginInfo;


// ****************************************************************************
//  Class: ViewerPlotFactory
//
//  Purpose:
//    ViewerPlotFactory is a factory for creating plots.  It also has
//    methods for manipulating the plot attributes.  When the class is
//    instantiated it registers all the different plot types.
//
//  Programmer: Eric Brugger
//  Creation:   August 23, 2000
//
//  Modifications:
//    Brad Whitlock, Fri Dec 8 17:16:55 PST 2000
//    I modified the CreatePlot method to accept more arguments.
//
//    Eric Brugger, Wed Dec 20 10:25:16 PST 2000
//    I modified the plot factory so that the constructor registers the
//    various plot types, instead of relying on the creator of the factory
//    to do so.  I also added methods for retrieving the default and client
//    attribute subjects for a given plot type.
//
//    Eric Brugger, Thu Mar  8 15:26:36 PST 2001
//    I modified the class to use the plot plugin manager.
//
//    Jeremy Meredith, Thu Jul 26 09:53:57 PDT 2001
//    Renamed plugin info to include the word "plot".
//
//    Jeremy Meredith, Fri Sep 28 13:47:32 PDT 2001
//    Removed the general plugin info since the viewer info is derived
//    from it now.
//
//    Brad Whitlock, Fri Apr 4 10:27:40 PDT 2003
//    I added nStates to the argument list.
//
// ****************************************************************************

class VIEWER_API ViewerPlotFactory
{
  public:
    ViewerPlotFactory();
    virtual ~ViewerPlotFactory();

    int GetNPlotTypes() const;

    ViewerPlot *CreatePlot(const int type, const char *hostName,
                           const char *databaseName, const char *var,
                           avtSILRestriction_p silr,
                           const int time0, const int time1,
                           const int nStates) const;

    AttributeSubject *GetDefaultAtts(const int type) const;
    AttributeSubject *GetClientAtts(const int type) const;

    void SetClientAttsFromDefault(const int type);
    void SetDefaultAttsFromClient(const int type);

  private:
    int                   nTypes;
    ViewerPlotPluginInfo  **viewerPluginInfo;
};

#endif
