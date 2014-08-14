/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
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

// ************************************************************************* //
//                               ViewerSubject.h                             //
// ************************************************************************* //

#ifndef VIEWER_SUBJECT_H
#define VIEWER_SUBJECT_H

#define VIEWER_NUM_ENGINE_RESTARTS        2

#include <viewer_exports.h>
#include <QObject>

#include <ViewerBase.h>
#include <ParentProcess.h>
#include <ViewerMasterXfer.h>
#include <VisWindowTypes.h>
#include <avtTypes.h>
#include <EngineKey.h>
#include <vectortypes.h>
#include <string>
#include <map>

class QApplication;
class QSocketNotifier;
class QTimer;

class BufferConnection;
class DataNode;
class HostProfileList;
class SILAttributes;
class ViewerActionBase;
class ViewerClientConnection;
class ViewerCommandFromSimObserver;
class ViewerConfigManager;
class ViewerMessageBuffer;
class ViewerMetaDataObserver;
class ViewerMethods;
class ViewerObserverToSignal;
class ViewerOperatorFactory;
class ViewerPlotFactory;
class ViewerSILAttsObserver;
class ViewerStateBuffered;
class ViewerWindow;
class avtDatabaseMetaData;
class avtDefaultPlotMetaData;
class SharedDaemon;

// ****************************************************************************
//  Class: ViewerSubject
//
//  Purpose:
//      ViewerSubject is the subject class for the Viewer proxy class.  It
//      implements all the methods in the Viewer class within the viewer
//      process.
//
//  Programmer: Eric Brugger
//  Creation:   August 9, 2000
//
//  Modifications:
//    Brad Whitlock, Thu Oct 5 19:38:33 PST 2000
//    I removed the SocketConnections since they are now in ParentProcess.
//
//    Brad Whitlock, Wed Nov 1 14:49:18 PST 2000
//    I made the class inherit from QObject so we can use Qt signals & slots.
//
//    Hank Childs, Fri Nov 10 13:19:59 PST 2000
//    Added TogglePerspective.
//
//    Eric Brugger, Fri Feb 23 12:27:00 PST 2001
//    I added the methods RecenterView and ToggleAutoCenterMode.
//
//    Brad Whitlock, Thu Apr 19 11:04:31 PDT 2001
//    Added methods to handle window iconification.
//
//    Brad Whitlock, Mon Apr 23 13:10:33 PST 2001
//    Added message-related things.
//
//    Brad Whitlock, Mon Apr 30 12:29:40 PDT 2001
//    Added methods to close and interrupt engines. Also added status methods.
//
//    Brad Whitlock, Wed Jun 13 17:33:10 PST 2001
//    Added a method to update color tables.
//
//    Kathleen Bonnell, Mon Jun 18 14:56:09 PDT 2001 
//    I added AnnotationAttributes, and method ApplyAnnotationAttributes. 
//
//    Jeremy Meredith, Tue Jul  3 15:14:19 PDT 2001
//    Added ReadFromParentAndCheckForInterruption, renamed ReadFromParent
//    to ReadFromParentAndProcess.
//
//    Brad Whitlock, Fri May 25 16:37:10 PST 2001
//    Made viewer quit if the gui dies.
//
//    Brad Whitlock, Thu Jun 21 13:50:02 PST 2001
//    Added a method to set SIL restrictions.
//
//    Jeremy Meredith, Fri Jul 20 11:24:21 PDT 2001
//    Added 'shift' to correct some window managers.
//
//    Brad Whitlock, Wed Jul 18 09:17:36 PDT 2001
//    Added methods to set the view.
//
//    Brad Whitlock, Tue Aug 14 17:32:21 PST 2001
//    Added methods to reset plot and operator attributes.
//
//    Eric Brugger, Tue Aug 21 10:13:26 PDT 2001
//    Removed the ViewCallback method.
//
//    Brad Whitlock, Thu Aug 30 08:59:08 PDT 2001
//    Removed the annotation attributes. Added methods to set default/reset
//    the annotation attributes.
//
//    Jeremy Meredith, Wed Sep  5 14:05:29 PDT 2001
//    Added plugin manager attributes.
//
//    Brad Whitlock, Wed Sep 5 00:29:10 PDT 2001
//    I added appearance attributes.
//
//    Jeremy Meredith, Fri Sep 14 13:48:47 PDT 2001
//    I added 'preshift' as a further correction for window managers.
//
//    Brad Whitlock, Fri Sep 21 13:24:04 PST 2001
//    I added a new Status method.
//
//    Brad Whitlock, Fri Sep 14 14:06:25 PST 2001
//    I added lighting support and changed the animation to use a different
//    style of timer.
//
//    Eric Brugger, Mon Oct 29 09:46:41 PST 2001
//    I removed the timer and work process.  I also added the method
//    ProcessFromParent.
//
//    Eric Brugger, Mon Nov 19 15:39:43 PST 2001
//    I added SetAnimationAttributes.
//
//    Brad Whitlock, Mon Sep 17 11:01:36 PDT 2001
//    I added a mechanism to send a synchronize event back to the client. I
//    added methods for suppressing redraws and for redrawing the window.
//
//    Kathleen Bonnell, Tue Nov 17 16:06:00 PST 200
//    I added ClearPickPoints and method sfor setting pick attributes.
//
//    Brad Whitlock, Fri Jan 4 17:26:58 PST 2002
//    I added a flag that determines whether or not interruption is enabled.
//
//    Brad Whitlock, Tue Jan 29 16:21:22 PST 2002
//    I added a new method to handle the SetWindowArea RPC.
//
//    Brad Whitlock, Wed Feb 20 14:30:35 PST 2002
//    I added a method for printing windows.
//
//    Brad Whitlock, Wed Mar 6 16:24:57 PST 2002
//    I added methods for handling file replacement and overlay. I also moved
//    all the infrastructure for default SIL restrictions to ViewerPlotList.
//
//    Jeremy Meredith, Tue Apr 16 11:25:54 PDT 2002
//    Added ability to process a few command line args before anything happens.
//    Added flag to disable reading of config files.
//
//    Jeremy Meredith, Tue Apr 23 15:51:24 PDT 2002
//    Added a nowin flag to disable window metric testing.
//
//    Kathleen Bonnell, Fri Apr 25 10:44:46 PST 2001
//    Add ClearRefLines.
//
//    Brad Whitlock, Mon May 6 17:32:48 PST 2002
//    Added some new methods that expose the popup menu's function to the
//    viewer client.
//
//    Jeremy Meredith, Wed May  8 15:39:06 PDT 2002
//    Added keyframe attributes.
//
//    Hank Childs, Wed May 29 08:42:29 PDT 2002
//    Added ToggleSpinMode.
//
//    Brad Whitlock, Thu Jun 27 16:21:02 PST 2002
//    Added methods to handle copying of window attributes.
//
//    Hank Childs, Mon Jul 15 09:45:42 PDT 2002
//    Allow for setting view based on different flavors of extents.
//
//    Brad Whitlock, Mon Jul 29 15:18:15 PST 2002
//    I added a method to reopen a database and a method to clear the
//    cache of a compute engine.
//
//    Brad Whitlock, Mon Aug 19 16:00:43 PST 2002
//    I removed the animationAtts member since it was not used.
//
//    Brad Whitlock, Thu Sep 19 13:30:02 PST 2002
//    I added methods to handle a couple new RPCs.
//
//    Brad Whitlock, Fri Sep 6 14:12:18 PST 2002
//    I added methods to handle queries.
//
//    Brad Whitlock, Fri Sep 27 15:29:39 PST 2002
//    I added support for a launch progress window.
//
//    Jeremy Meredith, Thu Oct 24 16:03:39 PDT 2002
//    Added material options.
//
//    Brad Whitlock, Mon Nov 11 11:49:27 PDT 2002
//    Added locking time and tools.
//
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002  
//    Removed SetDefault/Reset and Set for PickAttributes. 
//
//    Eric Brugger, Mon Nov 18 13:44:04 PST 2002
//    I added SetPlotFrameRange and DeletePlotKeyframe.
//
//    Eric Brugger, Mon Nov 25 09:41:36 PST 2002
//    I moved the keyframe attribute subject to the window manager.  I added
//    SetKeyframeAttributes.
//
//    Brad Whitlock, Mon Dec 2 13:36:18 PST 2002
//    I removed the ViewerColorTables member.
//
//    Jeremy Meredith, Thu Dec 19 12:08:46 PST 2002
//    Added support for launching engines from the command line.
//
//    Eric Brugger, Mon Dec 30 11:29:49 PST 2002
//    I added SetPlotDatabaseState and DeletePlotDatabaseKeyframe.
//
//    Eric Brugger, Fri Jan  3 16:00:27 PST 2003
//    I added ClearViewKeyframes, DeleteViewKeyframe, SetViewKeyframe and
//    ToggleCameraViewMode.
//
//    Brad Whitlock, Thu Jan 9 11:48:50 PDT 2003
//    I made the borders, shift, preshift, and geometry members be strings
//    instead of character arrays to avoid bad memory problems when closing
//    VisIt.
//
//    Brad Whitlock, Mon Jan 13 08:51:01 PDT 2003
//    I added OpenMDServer.
//
//    Eric Brugger, Tue Jan 28 14:23:23 PST 2003
//    I added MovePlotKeyframe, MovePlotDatabaseKeyframe and MoveViewKeyframe.
//
//    Brad Whitlock, Wed Feb 5 10:06:23 PDT 2003
//    I removed a lot of the slots since I totally changed how the popup
//    menu works when I added a toolbar to the windows.
//
//    Brad Whitlock, Thu Feb 27 11:19:47 PDT 2003
//    I added support for processing special opcodes, which are commands from
//    the client that are executed immediately.
//
//    Kathleen Bonnell, Wed Feb 19 11:56:32 PST 2003
//    Added SetGlobalLineoutAttributes.
//
//    Brad Whitlock, Wed Mar 12 10:32:35 PDT 2003
//    I made window iconification a special opcode.
//
//    Brad Whitlock, Thu Mar 20 13:07:08 PST 2003
//    I removed some methods related to plots and operators that I turned
//    into actions.
//
//    Eric Brugger, Fri Apr 18 12:43:54 PDT 2003 
//    I replaced ToggleAutoCenterMode with ToggleMaintainViewMode.
//
//    Brad Whitlock, Thu May 15 13:36:48 PST 2003
//    I added OpenDatabaseHelper.
//
//    Brad Whitlock, Fri May 16 14:52:32 PST 2003
//    I added configFileName.
//
//    Kathleen Bonnell, Tue Jul  1 09:21:57 PDT 2003 
//    Added SetPickAttributes.
//
//    Brad Whitlock, Mon Jun 30 12:20:37 PDT 2003
//    I added methods that allow more information to be saved to a config file.
//
//    Brad Whitlock, Tue Jul 1 16:59:00 PST 2003
//    Added ExportColorTable.
//
//    Brad Whitlock, Wed Jul 9 12:33:40 PDT 2003
//    Added ExportEntireState and ImportEntireState.
//
//    Eric Brugger, Wed Aug 20 11:09:33 PDT 2003
//    I added SetViewCurve.
//
//    Walter Herrera, Tue Set 9 10:31:04 PDT 2003
//    Added CreateAttributesDataNode
//
//    Jeremy Meredith, Fri Sep 26 12:18:58 PDT 2003
//    Added another flag for stereo so we could turn its default value
//    to 'true' if "-stereo" was on the command line.
//
//    Kathleen Bonnell, Wed Nov 26 14:35:29 PST 2003 
//    Added ResetPickAttributes. 
//    
//    Brad Whitlock, Wed Oct 29 11:00:55 PDT 2003
//    Added methods to handle several new annotation RPCs.
//
//    Kathleen Bonnell, Wed Dec 17 14:45:22 PST 2003 
//    Added SetDefaultPickAttributes, ResetPickLetter. 
//
//    Brad Whitlock, Thu Dec 18 12:25:14 PDT 2003
//    I added the processingFromParent member. I also changed some methods
//    that deal with processing input to be public slots.
//
//    Brad Whitlock, Thu Jan 29 22:20:56 PST 2004
//    I added ActivateDatabase, CheckForNewStates, CreateDatabaseCorrelation,
//    AlterDatabaseCorrelation, DeleteDatabaseCorrelation, and CloseDatabase.
//
//    Brad Whitlock, Thu Feb 26 13:32:43 PST 2004
//    Added ClearCacheForAllEngines.
//
//    Brad Whitlock, Fri Mar 12 12:08:59 PDT 2004
//    I added SendKeepAlives.
//
//    Jeremy Meredith, Tue Mar 23 14:36:37 PST 2004
//    Added engineParallelArguments.
//
//    Kathleen Bonnell, Wed Mar 24 07:13:33 PST 2004 
//    Added  QueryOverTimeAttributes methods.
//
//    Eric Brugger, Mon Mar 29 14:22:55 PST 2004
//    I added ToggleMaintainDataMode.
//
//    Brad Whitlock, Mon Aug 2 15:35:52 PST 2004
//    Added internal slot function EnableSocketSignals.
//
//    Kathleen Bonnell, Thu Aug  5 08:34:15 PDT 2004 
//    Added ResetLineoutColor.
//
//    Kathleen Bonnell, Wed Aug 18 09:28:51 PDT 2004 
//    Added InteractorAttributes methods.
//
//    Jeremy Meredith, Wed Aug 25 10:47:18 PDT 2004
//    Added ability to read information (both in general and specifically
//    for new metadata and SIL atts) from an engine.  This was needed for
//    simulations.
//
//    Brad Whitlock, Thu Feb 3 11:09:27 PDT 2005
//    Added int return value for OpenDatabaseHelper.
//
//    Mark C. Miller, Tue Mar  8 18:06:19 PST 2005
//    Added ProcessAttributees.
//
//    Brad Whitlock, Fri Mar 18 16:37:51 PST 2005
//    I removed ToggleLockTime.
//
//    Jeremy Meredith, Mon Apr  4 17:36:13 PDT 2005
//    Added SendSimulationCommand.
//
//    Brad Whitlock, Thu Apr 14 16:30:48 PST 2005
//    Added PostponeViewerRPC.
//
//    Hank Childs, Wed May 25 14:46:26 PDT 2005
//    Added UpdateDBPluginInfo.
//
//    Mark C. Miller, Tue May 31 20:12:42 PDT 2005
//    Added SetTryHarderCyclesTimes
//
//    Brad Whitlock, Tue May 3 15:14:16 PST 2005
//    Added new members that allow the viewer to talk to multiple clients. Part
//    of my new coding allowed me to remove the hack method that used to allow
//    the engine manager to check for interruption. Added movieAtts so they
//    can be kept in sync between multiple clients.
//
//    Hank Childs, Mon Jul 18 16:00:32 PDT 2005
//    Add qt_argv, which cannot be deleted until the program exits.
//
//    Kathleen Bonnell, Wed Jul 27 15:47:34 PDT 2005 
//    Added SuppressQueryOutput. 
//
//    Mark C. Miller, Wed Nov 16 10:46:36 PST 2005
//    Added mesh management attributes rpcs 
//
//    Brad Whitlock, Thu Nov 17 17:06:30 PST 2005
//    Added MoveWindow, ResizeWindow, MoveAndResizeWindow.
//
//    Hank Childs, Mon Feb 13 21:59:53 PST 2006
//    Added ConstructDDF.
//
//    Brad Whitlock, Thu May 11 15:05:17 PST 2006
//    Added ErrorClear.
//
//    Kathleen Bonnell, Tue Jun 20 16:02:38 PDT 2006 
//    Added UpatePlotInfoAtts. 
//
//    Jeremy Meredith, Mon Aug 28 16:55:01 EDT 2006
//    Added ability to force using a specific plugin when opening a file.
//
//    Hank Childs, Thu Jan 11 15:33:07 PST 2007
//    Added return value for OpenDatabase so that we can do a better job
//    when opening a bad file.
//
//    Brad Whitlock, Thu Jan 25 14:21:23 PST 2007
//    Added engineCommandObserver and the HandleCommandFromSimulation method.
//
//    Brad Whitlock, Tue Feb 13 14:00:29 PST 2007
//    I made it use ViewerState for all of its state objects.
//
//    Brad Whitlock, Fri Mar 9 16:25:33 PST 2007
//    Added HandleRequestMetaData.
//
//    Brad Whitlock, Mon Jun 11 11:56:37 PDT 2007
//    Added InterpretCommands method and some helper methods.
//
//    Cyrus Harrison, Tue Sep 18 11:13:27 PDT 2007
//    Added SetQueryFloatFormat
//
//    Kathleen Bonnell, Tue Oct  9 14:40:10 PDT 2007 
//    Added SetCreateMeshQualityExpressions and
//    SetCreateTimeDerivativeExpressions
//
//    Cyrus Harrison, Wed Nov 28 12:03:27 PST 2007
//    Added SetCreateVectorMagnitudeExpressions
//
//    Jeremy Meredith, Wed Jan 23 16:25:45 EST 2008
//    Added call to update the current default file opening options.
//
//    Brad Whitlock, Thu Jan 24 09:49:27 PST 2008
//    Added DeferCommandFromSimulation.
//
//    Brad Whitlock, Thu Jan 31 12:01:12 PST 2008
//    Added RemoveCrashRecoveryFile.
//
//    Jeremy Meredith, Mon Feb  4 13:29:43 EST 2008
//    Added SetViewAxisArray.
//
//    Brad Whitlock, Wed Feb 13 14:12:15 PST 2008
//    Added argument to SetFromNode.
//
//    Cyrus Harrison,Thu Feb 21 15:04:14 PST 2008
//    Added SetSuppressMessages()
//
//    Brad Whitlock, Thu Apr 10 10:09:31 PDT 2008
//    Added applicationLocale.
//
//    Brad Whitlock, Thu Aug 14 14:48:53 PDT 2008
//    Refactored the ViewerSubject so it no longer contains the QApplication or
//    the event loop so it can be embedded in other Qt applications.
//
//    Hank Childs, Wed Jan 28 14:53:32 PST 2009
//    Add support for named selections.
//
//    Brad Whitlock, Tue Apr 14 11:23:44 PDT 2009
//    I moved a bunch of members to the base class' ViewerProperties so we 
//    can access their values throughout the viewer.
//
//    Jeremy Meredith, Wed Feb  3 15:35:08 EST 2010
//    Removed maintain data; moved maintain view from Global settings
//    (Main window) to per-window Window Information (View window).
//
//    Jeremy Meredith, Wed Apr 21 13:19:16 EDT 2010
//    Save a copy of the host profiles we loaded from the system
//    installation directory, so we know if a user changed anything.
//
//    Hank Childs, Sat Aug 21 14:05:14 PDT 2010
//    Rename ConstructDDF to ConstructDataBinning.
//
//    Kathleen Biagas, Fri Jun 17 16:30:22 PDT 2011
//    Database, Line and PointQuery methods replaced by Query.
//
//    Kathleen Biagas, Fri Jul 15 11:34:11 PDT 2011
//    Added GetQueryParameters.
//
//    Marc Durant, Thu Jan 12 12:36:00 MST 2012
//    Added ToggleAllowPopup.
//
//    Jonathan Byrd (Allinea Software), Sun Dec 18 2011
//    Added methods for connecting/disconnecting the viewer with DDT,
//    and to instruct DDT to focus on a specific domain.
//
//    Kathleen Biagas, Wed Aug  7 12:59:21 PDT 2013
//    Add SetPrecisionType.
//
//    Cameron Christensen, Tuesday, June 10, 2014
//    Added SetBackendTypeRPC.
//
// ****************************************************************************

class VIEWER_API ViewerSubject : public ViewerBase
{
    Q_OBJECT
public:
    ViewerSubject();
    ~ViewerSubject();

    void ProcessCommandLine(int argc, char **argv);
    void Connect(int *argc, char ***argv);
    void Initialize();
    void RemoveCrashRecoveryFile() const;

    void SetNowinMode(bool);
    bool GetNowinMode() const;

    ViewerState   *GetViewerDelayedState();
    ViewerMethods *GetViewerDelayedMethods();
public:
    // Methods that should not be called outside of the "viewer".
    // ************************************************************************
    ViewerPlotFactory *GetPlotFactory() const;
    ViewerOperatorFactory *GetOperatorFactory() const;

    void MessageRendererThread(const char *message);

    static void ProcessEventsCB(void *cbData);
    static bool LaunchProgressCB(void *data, int stage);
    static void SpecialOpcodeCallback(int opcode, void *data);

    void StartLaunchProgress();
    void EndLaunchProgress();
    void BlockSocketSignals(bool);

    void CreateNode(DataNode *node, bool detailed);
    bool SetFromNode(DataNode *node, const std::map<std::string,std::string>&,
                    const std::string &);

    void PostponeAction(ViewerActionBase *);

    // Callback function for opening processes via engine.
    static void OpenWithEngine(const std::string &remoteHost, 
                               const stringVector &args, void *data);

    void AddNewViewerClientConnection(ViewerClientConnection* newClient);

public slots:
    void ProcessFromParent();
private:
    void CreateState();
    void ConnectXfer();
    void ConnectObjectsAndHandlers();
    void ConnectConfigManager();

    void InitializePluginManagers();
    void LoadPlotPlugins();
    void LoadOperatorPlugins();
    void InformClientOfPlugins() const;
    void ProcessConfigFileSettings();
    void AddInitialWindows();
    void LaunchEngineOnStartup();

    void ProcessEvents();
    void ReadConfigFiles(int argc, char **argv);
    void CustomizeAppearance();
    void InitializeWorkArea();
    int  OpenDatabaseHelper(const std::string &db, int timeState,
                            bool loadDefaultPlots, bool udpateWindowInfo,
                            const std::string &forcedFileType = "");

    bool HasInterpreter() const;
    void InterpretCommands(const std::string &commands);

    DataNode *CreateAttributesDataNode(const avtDefaultPlotMetaData *) const;

    // RPC handler methods.
    void HandleViewerRPCEx();
    void Close();
    void ConnectToMetaDataServer();
    void IconifyAllWindows();
    void DeIconifyAllWindows();
    void ShowAllWindows();
    void HideAllWindows();
    bool OpenDatabase();
    void CloseDatabase();
    void ActivateDatabase();
    void CheckForNewStates();
    void ReOpenDatabase();
    void ReplaceDatabase();
    void OverlayDatabase();
    void HandleRequestMetaData();
    void CreateDatabaseCorrelation();
    void AlterDatabaseCorrelation();
    void DeleteDatabaseCorrelation();
    void OpenComputeEngine();
    void CloseComputeEngine();
    void InterruptComputeEngine();
    void OpenMDServer();
    void ClearCache();
    void ClearCacheForAllEngines();
    void ExportDatabase();
    void ConstructDataBinning();
    void UpdateDBPluginInfo();

    void SendSimulationCommand();

    void ApplyNamedSelection();
    void CreateNamedSelection();
    void DeleteNamedSelection();
    void LoadNamedSelection();
    void SaveNamedSelection();
    void SetNamedSelectionAutoApply();
    void UpdateNamedSelection();
    void UpdateNamedSelection(const std::string &selName, bool updatePlots, bool allowCache);
    void InitializeNamedSelectionVariables();

    void SetDefaultPlotOptions();
    void ResetPlotOptions();
    void AddInitializedOperator();
    void SetDefaultOperatorOptions();
    void ResetOperatorOptions();

    void PrintWindow();
    void SaveWindow();
    void ExportWindow();
    void ExportHostProfile();
    void Export();
    void UpdateColorTable();
    void ExportColorTable();
    void WriteConfigFile();
    void ExportEntireState();
    void ImportEntireState();
    void ImportEntireStateWithDifferentSources();

    void SetAnnotationAttributes();
    void SetDefaultAnnotationAttributes();
    void ResetAnnotationAttributes();
    void AddAnnotationObject();
    void HideActiveAnnotationObjects();
    void DeleteActiveAnnotationObjects();
    void RaiseActiveAnnotationObjects();
    void LowerActiveAnnotationObjects();
    void SetAnnotationObjectOptions();
    void SetDefaultAnnotationObjectList();
    void ResetAnnotationObjectList();

    void SetKeyframeAttributes();
    void SetViewAxisArray();
    void SetViewCurve();
    void SetView2D();
    void SetView3D();
    void ClearViewKeyframes();
    void DeleteViewKeyframe();
    void MoveViewKeyframe();
    void SetViewKeyframe();
    void SetViewExtentsType();
    void SetAppearanceAttributes();
    void DisableRedraw();
    void RedrawWindow();
    void ProcessExpressions();
    void SetLightList();
    void SetDefaultLightList();
    void ResetLightList();
    void SetAnimationAttributes();
    void SetWindowArea();

    void CopyViewToWindow();
    void CopyLightingToWindow();
    void CopyAnnotationsToWindow();
    void CopyPlotsToWindow();
    void SetRenderingAttributes();
    void SuppressQueryOutput();
    void SetQueryFloatFormat();
    void Query();
    void GetQueryParameters();
    void SetMaterialAttributes();
    void SetDefaultMaterialAttributes();
    void ResetMaterialAttributes();
    void SetGlobalLineoutAttributes();
    void SetPickAttributes();
    void SetDefaultPickAttributes();
    void ResetPickAttributes();
    void ResetPickLetter();
    void RenamePickLabel();
    void SetQueryOverTimeAttributes();
    void SetDefaultQueryOverTimeAttributes();
    void ResetQueryOverTimeAttributes();
    void ResetLineoutColor();
    void SetMeshManagementAttributes();
    void SetDefaultMeshManagementAttributes();
    void ResetMeshManagementAttributes();

    void SetInteractorAttributes();
    void SetDefaultInteractorAttributes();
    void ResetInteractorAttributes();

    void GetProcessAttributes();

    void SetTryHarderCyclesTimes();
    void SetTreatAllDBsAsTimeVarying();

    void SetCreateMeshQualityExpressions();
    void SetCreateTimeDerivativeExpressions();
    void SetCreateVectorMagnitudeExpressions();
    void SetPrecisionType();
    void SetBackendType();

    void MoveWindow();
    void MoveAndResizeWindow();
    void ResizeWindow();

    void OpenClient();

    void SetDefaultFileOpenOptions();
    void SetSuppressMessages();
    void BroadcastAdvanced(AttributeSubject *subj);

    void DDTFocus();
    void DDTConnect();
        
signals:
    void scheduleHeavyInitialization();
private slots:
    void HeavyInitialization();

    void AddInputToXfer(ViewerClientConnection *, AttributeSubject *subj);
    void ProcessSpecialOpcodes(int opcode);
    void DisconnectClient(ViewerClientConnection *client);
    void DiscoverClientInformation();
    void CreateViewerDelayedState();

    void HandleViewerRPC();
    void HandlePostponedAction();
    void HandleSync();
    void HandleClientMethod();
    void HandleClientInformation();
    void HandleMetaDataUpdated(const std::string &host, const std::string &db,
                               const avtDatabaseMetaData *md);
    void HandleSILAttsUpdated(const std::string &host, const std::string &db,
                              const SILAttributes *md);
    void DeferCommandFromSimulation(const EngineKey &key, const std::string &db,
                                    const std::string &command);
    void HandleCommandFromSimulation(const EngineKey &key, const std::string &db,
                                     const std::string &command);
    void HandleColorTable();
    void ProcessRendererMessage();
    void ReadFromParentAndProcess(int);
    void DelayedProcessSettings();
    void SendKeepAlives();
    void EnableSocketSignals();

    void ReadFromSimulationAndProcess(int);

    void ConnectWindow(ViewerWindow *win);

    void ToggleMaintainViewMode(int windowIndex = -1);
    void ToggleCameraViewMode(int windowIndex = -1);
    void ToggleLockTools(int windowIndex = -1);
    void ToggleAllowPopup(int windowIndex = -1);

    void CopyViewToWindow(int from, int to);
    void CopyLightingToWindow(int from, int to);
    void CopyAnnotationsToWindow(int from, int to);
    void CopyPlotsToWindow(int from, int to);

    void OpenDatabaseOnStartup();
    void OpenScriptOnStartup();
private:
    typedef std::vector<ViewerClientConnection *> ViewerClientConnectionVector;
    static void BroadcastToAllClients(void *, Subject *);

    QSocketNotifier       *checkParent;
    QTimer                *keepAliveTimer;
    bool                   launchingComponent;
    bool                   deferHeavyInitialization;
    bool                   heavyInitializationDone;
    bool                   interruptionEnabled;
    std::string            launchEngineAtStartup;
    std::string            openDatabaseOnStartup;
    std::string            openScriptOnStartup;
    bool                   blockSocketSignals;
    bool                   processingFromParent;
    int                    animationStopOpcode;
    int                    iconifyOpcode;
    std::string            interpretCommands;

    ViewerMasterXfer       xfer;
    ParentProcess         *parent;
    ViewerClientConnectionVector clients;
    BufferConnection      *inputConnection;
    ViewerStateBuffered   *viewerDelayedState;
    ViewerMethods         *viewerDelayedMethods;

    ViewerObserverToSignal *viewerRPCObserver;
    ViewerObserverToSignal *postponedActionObserver;
    ViewerObserverToSignal *syncObserver;
    ViewerObserverToSignal *clientMethodObserver;
    ViewerObserverToSignal *clientInformationObserver;
    ViewerObserverToSignal *colorTableObserver;

    ViewerPlotFactory     *plotFactory;
    ViewerOperatorFactory *operatorFactory;
    ViewerConfigManager   *configMgr;
    DataNode              *systemSettings;
    DataNode              *localSettings;

    HostProfileList       *originalSystemHostProfileList;

    ViewerMessageBuffer   *messageBuffer;

    std::map<int,EngineKey>                     simulationSocketToKey;
    std::map<EngineKey,QSocketNotifier*>        engineKeyToNotifier;
    std::map<EngineKey,ViewerMetaDataObserver*> engineMetaDataObserver;
    std::map<EngineKey,ViewerSILAttsObserver*>  engineSILAttsObserver;
    std::map<EngineKey,ViewerCommandFromSimObserver*>  engineCommandObserver;

    // note: we may only want to use the engineParallelArguments for
    // the first launch of an engine
    std::vector<std::string> engineParallelArguments;
    std::vector<std::string> unknownArguments;
    std::vector<std::string> clientArguments;
    SharedDaemon             *shared_viewer_daemon;
    size_t                   clientIds;
};

#endif
