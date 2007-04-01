// ************************************************************************* //
//                            ViewerEngineManager.h                          //
// ************************************************************************* //

#ifndef VIEWER_ENGINE_MANAGER_H
#define VIEWER_ENGINE_MANAGER_H
#include <viewer_exports.h>
#include <ViewerServerManager.h>
#include <SimpleObserver.h>
#include <avtDataObjectReader.h>
#include <vectortypes.h>
#include <EngineKey.h>
#include <map>

// Forward declarations.
class AttributeSubject;
class EngineList;
class EngineProxy;
class MaterialAttributes;
class PickAttributes;
class QueryAttributes;
class QueryOverTimeAttributes;
class ViewerPlot;
class AnnotationAttributes;
class WindowAttributes;

// ****************************************************************************
//  Class: ViewerEngineManager
//
//  Purpose:
//      ViewerEngineManager is a class that manages engines.  It has methods
//      for creating and closing engines.  It also has a method to get the
//      vtkDataSet for a plot.  It is modeled after the singleton class which
//      means that applications can't directly instantiate the class but can
//      only get a reference to the one instance through the Instance method.
//
//  Programmer: Eric Brugger
//  Creation:   September 22, 2000
//
//  Modifications:
//    Hank Childs, Thu Sep 28 22:16:55 PDT 2000
//    Made the GetDataSet routine return an avtDataset.
//
//    Hank Childs, Tue Oct 17 09:31:02 PDT 2000
//    Changed GetDataSet to GetDataSetReader.
//
//    Brad Whitlock, Mon Nov 27 17:26:20 PST 2000
//    Added code to set the debug level that gets passed to the engines.
//
//    Hank Childs, Fri Dec 29 08:39:27 PST 2000
//    Changed DataSetReader to DataObjectReader.
//
//    Jeremy Meredith, Fri Apr 20 10:35:47 PDT 2001
//    Added other arguments to engine creation.
//
//    Brad Whitlock, Tue Apr 24 12:52:53 PDT 2001
//    Added the RemoveEngine, GetEngineIndex methods.
//
//    Brad Whitlock, Mon Apr 30 12:35:42 PDT 2001
//    Added CloseEngine, InterruptEngine methods. Made the class inherit from
//    SimpleObserver so it can see status updates reported by the engine
//    proxies.
//
//    Jeremy Meredith, Tue Jul  3 15:13:47 PDT 2001
//    Added InterruptEngine(int) in case we know the index but no hostname.
//
//    Brad Whitlock, Mon Sep 24 13:59:44 PST 2001
//    Added a method to set the localHost variable.
//
//    Jeremy Meredith, Fri Nov  9 10:09:55 PST 2001
//    Added the ability to get an engine proxy by hostname.
//
//    Hank Childs, Tue Nov 20 14:51:59 PST 2001
//    Add a method to get an image for an existing plot.
//
//    Brad Whitlock, Fri May 3 16:07:41 PST 2002
//    Made it inherit from ViewerServerManager.
//
//    Brad Whitlock, Tue May 7 17:23:34 PST 2002
//    Added a couple protected methods.
//
//    Brad Whitlock, Tue Jul 30 13:20:27 PST 2002
//    I added a method to clear the cache for an engine.
//
//    Kathleen Bonnell, Wed Sep 18 10:33:31 PDT 2002 
//    I added a method to Query an engine and ReleaseData from an engine.
//
//    Jeremy Meredith, Thu Oct 24 16:15:11 PDT 2002
//    Added material options.
//
//    Brad Whitlock, Tue Dec 10 15:23:53 PST 2002
//    I added a method to tell the engine to open a database.
//
//    Eric Brugger, Thu Dec 19 11:00:07 PST 2002
//    I added keyframing support.
//
//    Jeremy Meredith, Thu Dec 19 12:17:52 PST 2002
//    Added ability to skip engine profile chooser window.  It is still
//    enabled by default.
//
//    Brad Whitlock, Fri Dec 27 11:59:08 PDT 2002
//    I added a vector of strings to hold the arguments that are used when
//    restarting a crashed engine.
//
//    Brad Whitlock, Tue Mar 25 14:12:34 PST 2003
//    I added the DefineVirtualDatabase rpc that allows us to define visit
//    file style databases without having to have a visit file on disk.
//
//    Mark C. Miller, 08Apr03
//    Added another GetDataObjectReader method for getting scalable rendered 
//    image from engine
//
//    Hank Childs, Thu Oct  2 16:18:11 PDT 2003
//    Allow queries to have multiple network ids.
//
//    Mark C. Miller, Wed Oct 29 16:56:14 PST 2003
//    Added method UseDataObjectReader
//
//    Mark C. Miller, Sat Jan 17 12:40:16 PST 2004
//    Removed 'numRestarts' data member. Made numRestarts argument to
//    CreateEngine default to -1 (e.g. unknown)
//
//    Brad Whitlock, Sat Jan 31 23:22:26 PST 2004
//    I removed a frame argument from GetDataObjectReader.
//
//    Brad Whitlock, Thu Feb 26 13:00:16 PST 2004
//    I added ClearCacheForMultipleEngines.
//
//    Hank Childs, Fri Mar  5 11:12:20 PST 2004
//    Add a format when opening a database.
//
//    Brad Whitlock, Fri Mar 12 11:39:27 PDT 2004
//    Added SendKeepAlives.
//
//    Jeremy Meredith, Mon Mar 22 09:20:58 PST 2004
//    Made CreateEngine return a success flag.
//
//    Jeremy Meredith, Thu Mar 25 15:56:48 PST 2004
//    Use a map of engines based on a key, and be aware of simulations.
//
//    Mark C. Miller, Mon Mar 29 14:52:08 PST 2004
//    Added new bool args to control annotations in SR mode on engine
//    GetDataObject reader used in external rendering and ExternalRender
//
//    Eric Brugger, Tue Mar 30 13:25:30 PST 2004
//    Added the plot data extents to MakePlot.
//
//    Kathleen Bonnell, Wed Mar 31 16:40:15 PST 2004 
//    Added CloneNetwork. 
//
//    Jeremy Meredith, Fri Apr  2 14:29:25 PST 2004
//    Made restartArguments be saved on a per-host (per-enginekey) basis.
//
//    Mark C. Miller, Wed Apr 14 16:41:32 PDT 2004
//    Added argument for extents type string to ExternalRender and 
//    SetWinAnnotAtts
//
// ****************************************************************************

class VIEWER_API ViewerEngineManager : public ViewerServerManager,
    public SimpleObserver
{
  public:
    virtual ~ViewerEngineManager();
    virtual void Update(Subject *TheChangedSubject);

    static ViewerEngineManager *Instance();

    bool CreateEngine(const EngineKey &ek,
                      const std::vector<std::string> &arguments,
                      bool  skipChooser=false,
                      int numRestarts=-1);
    bool ConnectSim(const EngineKey &ek,
                    const std::vector<std::string> &arguments,
                    const std::string &simHost, int simPort);
    void CloseEngines();
    void CloseEngine(const EngineKey &ek);
    void InterruptEngine(const EngineKey &ek);
    bool InExecute() const;
    void SendKeepAlives();

    bool ExternalRender(std::vector<const char*> pluginIDsList,
                        std::vector<EngineKey> engineKeysList,
                        intVector plotIdsList,
                        std::vector<const AttributeSubject *> attsList,
                        WindowAttributes winAtts,
                        AnnotationAttributes annotAtts,
                        std::string extStr,
                        bool& shouldTurnOfScalableRendering,
                        bool doAllAnnotations,
                        std::vector<avtImage_p>& imgList);

    avtDataObjectReader_p GetDataObjectReader(ViewerPlot *const plot);
    avtDataObjectReader_p UseDataObjectReader(ViewerPlot *const plot,
                                              bool turningOffScalableRendering);
    avtDataObjectReader_p GetDataObjectReader(bool sendZBuffer,
                                              const EngineKey &ek,
                                              intVector ids,
                                              bool doAllAnnotations);
    void GetImage(int index, avtDataObject_p &);
    void UpdatePlotAttributes(const std::string &,int index,AttributeSubject*);
    void ClearCacheForAllEngines();

    static EngineList *GetEngineList();

    // Global engine computation attributes
    static MaterialAttributes *GetMaterialClientAtts();
    static MaterialAttributes *GetMaterialDefaultAtts();
    static void SetClientMaterialAttsFromDefault();
    static void SetDefaultMaterialAttsFromClient();

    //
    // Engine RPCs
    //
    bool OpenDatabase( const EngineKey &ek, const char *format, 
                      const char *filename, int time=0);
    bool DefineVirtualDatabase(const EngineKey &ek, const char *format,
                               const char *dbName, const char *path, 
                               const stringVector &files, int time=0);
    bool ApplyOperator(const EngineKey &ek, const char *name,
                       const AttributeSubject *atts);
    bool MakePlot(const EngineKey &ek, const char *name,
                  const AttributeSubject *atts, const vector<double> &,
                  int *networkId);
    bool UseNetwork(const EngineKey &ek, int id);
    bool UpdatePlotAttributes(const EngineKey &ek, const char *name,
                              int id, const AttributeSubject *atts);
    bool Pick(const EngineKey &ek, const int nid,
              const PickAttributes *atts, PickAttributes &retAtts);
    bool StartPick(const EngineKey &ek,
                   const bool flag, const int nid);
    bool SetWinAnnotAtts(const EngineKey &ek, const WindowAttributes *wa,
                         const AnnotationAttributes *aa,
                         const std::string extStr);
    bool ClearCache(const EngineKey &ek, const char *dbName = 0);
    bool Query(const EngineKey &ek, const std::vector<int> &networkIds, 
               const QueryAttributes *atts, QueryAttributes &retAtts);
    bool ReleaseData(const EngineKey &ek, int id);
    bool CloneNetwork(const EngineKey &ek, int id, 
                      const QueryOverTimeAttributes *qatts);

  protected:
    ViewerEngineManager();
    void RemoveEngine(const EngineKey &ek, bool close);
    void RemoveFailedEngine(const EngineKey &ek);
    bool EngineExists(const EngineKey &ek) const;
    void UpdateEngineList();
    EngineProxy *GetEngine(const EngineKey &ek);
    void LaunchMessage(const EngineKey &ek) const;
    void BeginEngineExecute();
    void EndEngineExecute();

  private:
    typedef std::map<EngineKey,EngineProxy*> EngineMap;

    static ViewerEngineManager *instance;
    static EngineList          *clientEngineAtts;

    bool                       executing;
    int                        nEngines;
    EngineMap                  engines;

    std::map<EngineKey,stringVector>  restartArguments;

    // Global engine computation attributes
    static MaterialAttributes *materialClientAtts;
    static MaterialAttributes *materialDefaultAtts;
};

#endif
