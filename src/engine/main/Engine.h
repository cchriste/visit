/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef ENGINE_H
#define ENGINE_H

// RPCs
class ApplyOperatorRPC;
class ClearCacheRPC;
class CloneNetworkRPC;
class ConstructDDFRPC;
class DefineVirtualDatabaseRPC;
class ExecuteRPC;
class ExportDatabaseRPC;
class KeepAliveRPC;
class MakePlotRPC;
class OpenDatabaseRPC;
class PickRPC;
class QueryRPC;
class QuitRPC;
class ReadRPC;
class ReleaseDataRPC;
class RenderRPC;
class SetWinAnnotAttsRPC;
class SimulationCommand;
class StartPickRPC;
class StartQueryRPC;
class UpdatePlotAttsRPC;
class UseNetworkRPC;
class ProcInfoRPC;
class SimulationCommandRPC;
class SetEFileOpenOptionsRPC;

// Other classes
class avtDatabaseMetaData;
class Connection;
class LoadBalancer;
class MPIXfer;
class NetworkManager;
class NonBlockingRPC;
class Observer;
class ProcessAttributes;
class Xfer;

#include <vector>
#include <avtDataObjectWriter.h>
#include <ParentProcess.h>
#include <BufferConnection.h>

#ifdef PARALLEL
#include <parallel.h>
#endif

// ****************************************************************************
//  Class:  Engine
//
//  Purpose:
//    VisIt Computation Engine
//
//  Programmer:  Jeremy Meredith
//  Creation:    July  9, 2003
//
//  Modifications:
//    Jeremy Meredith, Wed Jan 14 10:38:48 PST 2004
//    Made some changes so that the engine can be used from within an
//    external event loop.
//
//    Brad Whitlock, Fri Mar 12 11:08:21 PDT 2004
//    Added KeepAliveRPC and SendKeepAliveReply.
//
//    Kathleen Bonnell, Wed Mar 31 16:53:03 PST 2004
//    Added CloneNetworkRPC.
//
//    Mark C. Miller, Mon May 24 18:36:13 PDT 2004
//    Added args to WriteData method to support checking if SR threshold is
//    exceeded
//
//    Hank Childs, Tue Jun  1 14:02:44 PDT 2004
//    Added Finalize method.
//
//    Mark C. Miller, Wed Jul  7 11:42:09 PDT 2004
//    Added NewHandler
//
//    Mark C. Miller, Wed Aug 11 23:42:18 PDT 2004
//    Added argument for cellCountMultiplier to WriteData
//
//    Jeremy Meredith, Wed Aug 25 11:58:44 PDT 2004
//    Added methods needed for Simulation support.
//
//    Jeremy Meredith, Mon Nov  1 17:17:59 PST 2004
//    Added a parallel ProcessInput for parallel Simulation support.  Made
//    some other small changes necessary for that to happen.
//
//    Hank Childs, Mon Feb 28 17:03:06 PST 2005
//    Added StartQueryRPC.
//
//    Mark C. Miller, Tue Mar  8 18:06:19 PST 2005
//    Added GetProcessAttributes
//
//    Hank Childs, Tue Mar 22 09:37:58 PST 2005
//    Remove unused data members that date back to old expression calls.
//
//    Jeremy Meredith, Mon Apr  4 15:59:49 PDT 2005
//    Added methods to control a simulation from VisIt.
//
//    Jeremy Meredith, Wed May 18 16:50:59 PDT 2005
//    Allow disabling of signal handlers.
//
//    Hank Childs, Thu May 26 11:50:45 PDT 2005
//    Add ExportDatabase.
//
//    Mark C. Miller, Wed Dec 14 17:19:38 PST 2005
//    Added compression bool to WriteData
//
//    Hank Childs, Thu Jan  5 14:24:06 PST 2006
//    Added a Boolean for whether or not we are using "-dump".
//
//    Hank Childs, Mon Feb 13 22:25:04 PST 2006
//    Added ConstructDDF.
//
//    Hank Childs, Tue Sep  5 10:45:13 PDT 2006
//    Make ResetTimeout be a public function, since query over time and
//    some of the line scan algorithms may take more than ten minutes.
//
//    Brad Whitlock, Thu Jan 25 15:26:45 PST 2007
//    Added SimulationInitiateCommand.
//
//    Jeremy Meredith, Wed Jan 23 16:50:36 EST 2008
//    Added setEFileOpenOptionsRPC.
//
//    Cyrus Harrison, Tue Feb 19 08:42:51 PST 2008
//    Removed shouldDoDashDump (flag now contained in avtDebugDumpOptions)
//
//    Brad Whitlock, Tue Jun 24 15:20:05 PDT 2008
//    Made plugin managers be members of the class instead of singletons.
//
//    Tom Fogal, Fri Jul 11 12:00:29 EDT 2008
//    Added a flag to determine whether or not we should use the IceT network
//    manager.  This should be set if we're given the `-icet' command line
//    parameter.
//
//    Tom Fogal, Tue Aug  5 14:26:03 EDT 2008
//    Added a string to store arguments to the X server.
//
// ****************************************************************************

class Engine
{
  public:
                   ~Engine();
    static Engine  *Instance();

    // Initialization routines
    void            Initialize(int *argc, char **argv[], bool sigs);
    void            Finalize(void);
    bool            ConnectViewer(int *argc, char **argv[]);
    void            SetUpViewerInterface(int *argc, char **argv[]);

    // Needed for simulations
    void            PopulateSimulationMetaData(const std::string &db,
                                               const std::string &fmt);
    void            SimulationTimeStepChanged();
    void            SimulationInitiateCommand(const char *);
    void            SetSimulationCommandCallback(void(*)(const char*,
                                                       int,float,const char*));
    void            ExecuteSimulationCommand(const std::string&,
                                             int,float,const std::string&);
    static void     DisconnectSimulation();

    // Two event loops
    bool            EventLoop();
    void            PAR_EventLoop();

    // Get the network manager
    NetworkManager *GetNetMgr()                   { return netmgr; }

    // Methods needed for an external event loop
    int             GetInputSocket();
    void            ProcessInput();
    void            PAR_ProcessInput();

    // Method to write data back to the viewer
    void            WriteData(NonBlockingRPC *, avtDataObjectWriter_p &,
                        bool useCompression=false,
                        bool respondWithNull=false, int scalableThresold=-1,
                        bool *scalableThresholdExceeded=0,
                        int currentTotalGlobalCellCount=0,
                        float cellCountMultiplier=1.0,
                        int* currentNetworkGlobalCellCount=0);
    void            SendKeepAliveReply();

    // Tell the engine whether or not fatal exceptions have occurred
    void            SetNoFatalExceptions(bool nf) {noFatalExceptions = nf;}

    // Method to get unix process information
    ProcessAttributes *GetProcessAttributes();

    // Various callbacks
    static bool     EngineAbortCallback(void *);
    static bool     EngineAbortCallbackParallel(void*, bool);
    static void     EngineInitializeProgressCallback(void *, int);
    static void     EngineUpdateProgressCallback(void*, const char*,
                                                 const char*, int,int);
    static void     EngineWarningCallback(void *, const char *);
    void            ResetTimeout(int timeout);
    void            SetOverrideTimeout(int mins)
                    {
                        overrideTimeoutEnabled = true;
                        overrideTimeoutMins = mins; 
                    }
    bool            IsIdleTimeoutEnabled(void) { return idleTimeoutEnabled; }

    // Internal methods
  protected:
                    Engine();
    void            ProcessCommandLine(int argc, char *argv[]);

    static void     AlarmHandler(int signal);
    static void     NewHandler(void);

  protected:
    // The singleton object
    static Engine     *instance;

    // The Viewer
    ParentProcess      theViewer;

    // The destination machine's type representation.
    TypeRepresentation destinationFormat;

    // Here's the network manager!
    NetworkManager    *netmgr;
    std::string        pluginDir;

    // The connection
    Connection        *vtkConnection;

    // Flag for whether there were any fatal exceptions
    bool               noFatalExceptions;

    // Current timeouts in minutes.
    int                idleTimeoutMins;
    int                executionTimeoutMins;
    bool               idleTimeoutEnabled;
    int                overrideTimeoutMins;
    bool               overrideTimeoutEnabled;

    // The load balancer
    LoadBalancer       *lb;


    // The Xfer object
#ifdef PARALLEL
    MPIXfer            *xfer;
    PAR_StateBuffer     par_buf;
    BufferConnection    par_conn;
#else
    Xfer               *xfer;
#endif

    // The list of all RPC executors (needed to clean up memory)
    std::vector<Observer*> rpcExecutors;

    // The RPCs
    QuitRPC                  *quitRPC;
    KeepAliveRPC             *keepAliveRPC;
    ReadRPC                  *readRPC;
    ApplyOperatorRPC         *applyOperatorRPC;
    MakePlotRPC              *makePlotRPC;
    UseNetworkRPC            *useNetworkRPC;
    UpdatePlotAttsRPC        *updatePlotAttsRPC;
    PickRPC                  *pickRPC;
    StartPickRPC             *startPickRPC;
    StartQueryRPC            *startQueryRPC;
    ExecuteRPC               *executeRPC;
    ClearCacheRPC            *clearCacheRPC;
    QueryRPC                 *queryRPC;
    ReleaseDataRPC           *releaseDataRPC;
    OpenDatabaseRPC          *openDatabaseRPC;
    DefineVirtualDatabaseRPC *defineVirtualDatabaseRPC;
    RenderRPC                *renderRPC;
    SetWinAnnotAttsRPC       *setWinAnnotAttsRPC;
    CloneNetworkRPC          *cloneNetworkRPC;
    ProcInfoRPC              *procInfoRPC;
    SimulationCommandRPC     *simulationCommandRPC;
    ExportDatabaseRPC        *exportDatabaseRPC;
    ConstructDDFRPC          *constructDDFRPC;
    SetEFileOpenOptionsRPC   *setEFileOpenOptionsRPC;

    // The metadata, filename, format, control data for a simulation
    std::string               filename;
    std::string               format;
    avtDatabaseMetaData      *metaData;
    SILAttributes            *silAtts;
    SimulationCommand        *commandFromSim;
    void                    (*simulationCommandCallback)(const char*,
                                                        int,float,const char*);

    // unix process attributes
    ProcessAttributes        *procAtts;

 private:
    bool                      useIceT;
    std::string               X_Args;
};

#endif
