#ifndef NETWORK_MANAGER_H
#define NETWORK_MANAGER_H

#include <avtDataObjectWriter.h>
#include <avtPlot.h>
#include <AnnotationAttributes.h>
#include <WindowAttributes.h>
#include <vectortypes.h>
#include <string>
#include <deque>
#include <map>

class AttributeGroup;
class CompactSILRestrictionAttributes;
class LoadBalancer;
class DataNetwork;
class Netnode;
class NetnodeDB;
class PickAttributes;
class QueryAttributes;
class MaterialAttributes;
class VisWindow;

// ****************************************************************************
//  Class: NetworkManager
//
//  Purpose:
//      Handles creation and caching of AVT networks.
//
//  Programmer: Jeremy Meredith
//  Creation:   September 29, 2000
//
//  Modifications:
//
//    Kathleen Bonnell, Thu Oct 12 12:50:27 PDT 2000
//    Added 'AddOnionPeel'.
//
//    Kathleen Bonnell, Tue Nov  7 15:17:38 PST 2000
//    Added 'AddRelevantPoints'.
//
//    Kathleen Bonnell, Fri Nov 17 16:33:40 PST 2000 
//    Added 'AddFilledBoundary'.
//
//    Jeremy Meredith, Tue Dec 12 13:52:48 PST 2000
//    Added 'AddMaterialSelect'.
//
//    Kathleen Bonnell, Wed Feb 28 15:04:30 PST 2001 
//    Added 'MakeContourPlot'.
//
//    Jeremy Meredith, Thu Mar  1 13:44:46 PST 2001
//    Removed all the AddXXXFilter methods and replaced them with a single
//    AddFilter factory-ish method.
//
//    Jeremy Meredith, Fri Mar  2 13:04:23 PST 2001
//    Replaced all the MakeXXXPlot methods with a single MakePlot
//    factory-ish method.
//
//    Jeremy Meredith, Sun Mar  4 16:59:11 PST 2001
//    Made AddFilter and MakePlot take a string instead of an enumerated
//    type.  This is to use the new PluginManager.
//
//    Jeremy Meredith, Fri Nov  9 10:23:13 PST 2001
//    Added EndNetwork, UseNetwork, and SetWindowAttributes.
//    Added a netMRU since the netPool is persistent.
//
//    Hank Childs, Fri Nov 30 12:37:26 PST 2001
//    Added UpdatePlotAtts.
//
//    Kathleen Bonnell, Mon Nov 19 13:17:52 PST 2001 
//    Added methods StartPickMode, StopPickMode, Pick, and member
//    requireOriginalCells.
//    
//    Hank Childs, Fri Dec 14 17:39:36 PST 2001
//    Added support for a more compact for of sil restrictions.
//
//    Hank Childs, Mon Jan  7 18:14:08 PST 2002
//    Clean up memory better.
//
//    Sean Ahern, Fri Apr 19 13:45:32 PDT 2002
//    Removed ApplyUnaryOperator.  Added ApplyNamedFunction.
//
//    Brad Whitlock, Fri Jun 28 14:02:38 PST 2002
//    Renamed Network class to DataNetwork.
//
//    Hank Childs, Wed Sep 11 09:32:21 PDT 2002
//    Add some support for cleaning up old networks.
//
//    Kathleen Bonnell, Mon Sep 16 14:28:09 PDT 2002   
//    Add Query. 
//
//    Jeremy Meredith, Thu Oct 24 16:03:39 PDT 2002
//    Added material options.
//
//    Mark C. MIller, Mon Dec  9 17:19:02 PST 2002
//    Added VisWindow data member
//
//    Brad Whitlock, Tue Dec 10 14:52:05 PST 2002
//    Overloaded AddDB.
//
//    Sean Ahern, Mon Dec 23 11:23:13 PST 2002
//    Rearranged the AddDB code into GetDBFromCache to reduce code
//    duplication.  Renamed AddDB to StartNetwork.
//
//    Brad Whitlock, Tue Mar 25 13:43:10 PST 2003
//    Added DefineDB.
//
//    Mark C. Miller, 15Jul03
//    Added method to set annotation attributes
//
//    Jeremy Meredith, Mon Sep 15 17:15:09 PDT 2003
//    Removed SetFinalVariableName.
//
//    Hank Childs, Thu Oct  2 16:31:08 PDT 2003
//    Allow queries to involve multiple networks.
//
// ****************************************************************************
class NetworkManager
{
   typedef std::map<std::string, stringVector> StringVectorMap;
 public:
                  NetworkManager(void);
                 ~NetworkManager(void);
    void          ClearAllNetworks(void);

    NetnodeDB*    GetDBFromCache(const string &filename, int time);
    void          StartNetwork(const std::string&, const std::string &,
                               int,
                               const CompactSILRestrictionAttributes &,
                               const MaterialAttributes &);
    void          DefineDB(const std::string &, const std::string &,
                           const stringVector &, int);
    void          AddFilter(const std::string&,
                            const AttributeGroup* = NULL,
                            const unsigned int ninputs = 1);
    void          MakePlot(const std::string&, const AttributeGroup* = NULL);
    int           EndNetwork(void);

    void          UseNetwork(int);
    avtPlot_p     GetPlot(void);
    int           GetCurrentNetworkId(void);
    int           GetTotalGlobalCellCounts(void) const;
    void          DoneWithNetwork(int);

    void          UpdatePlotAtts(int, const AttributeGroup *);

    void          SetWindowAttributes(const WindowAttributes&);
    void          SetAnnotationAttributes(const AnnotationAttributes&);

    void          SetLoadBalancer(LoadBalancer *lb) {loadBalancer = lb;};

    avtDataObjectWriter_p GetOutput(bool respondWithNullData,
                                    bool calledForRender);
    avtDataObjectWriter_p Render(intVector networkIds, bool getZBuffer);
 
    void          StartPickMode(void);
    void          StopPickMode(void);

    void          Pick(const int, PickAttributes *);
    void          Query(const std::vector<int> &, QueryAttributes*);

 private:
    std::vector<DataNetwork*>   networkCache;
    std::vector<int>            globalCellCounts;
    std::deque<DataNetwork*>    networkMRU;
    std::vector<NetnodeDB*>     databaseCache;
    StringVectorMap             virtualDatabases;

    std::vector<Netnode*>       workingNetnodeList;
    DataNetwork                *workingNet;
    std::vector<std::string>    nameStack;

    bool                        requireOriginalCells;
    LoadBalancer               *loadBalancer;
    WindowAttributes            windowAttributes;
    AnnotationAttributes        annotationAttributes;
    VisWindow                  *viswin;
};

#endif
