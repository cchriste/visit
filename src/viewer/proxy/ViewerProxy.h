/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                                ViewerProxy.h                              //
// ************************************************************************* //

#ifndef VIEWER_PROXY_H
#define VIEWER_PROXY_H
#include <viewerproxy_exports.h>
#include <SimpleObserver.h>
#include <avtSILRestriction.h>
#include <vectortypes.h>

class avtDatabaseMetaData;
class AnimationAttributes;
class AnnotationAttributes;
class AnnotationObjectList;
class AppearanceAttributes;
class AttributeSubject;
class ClientMethod;
class ClientInformation;
class ClientInformationList;
class ColorTableAttributes;
class Connection;
class ConstructDDFAttributes;
class DatabaseCorrelationList;
class DBPluginInfoAttributes;
class ExportDBAttributes;
class ExpressionList;
class EngineList;
class GlobalAttributes;
class GlobalLineoutAttributes;
class HostProfileList;
class InteractorAttributes;
class InternalSILObserver;
class KeyframeAttributes;
class LightList;
class MaterialAttributes;
class MeshManagementAttributes;
class MessageAttributes;
class MovieAttributes;
class ParentProcess;
class PickAttributes;
class PlotInfoAttributes;
class PlotList;
class PluginManagerAttributes;
class PostponedAction;
class PrinterAttributes;
class ProcessAttributes;
class QueryAttributes;
class QueryOverTimeAttributes;
class QueryList;
class RenderingAttributes;
class RemoteProcess;
class SaveWindowAttributes;
class SILAttributes;
class StatusAttributes;
class SyncAttributes;
class ViewCurveAttributes;
class View2DAttributes;
class View3DAttributes;
class ViewerRPC;
class WindowInformation;
class Xfer;

// ****************************************************************************
//  Class: ViewerProxy
//
//  Purpose:
//    ViewerProxy is a proxy class for creating and controlling a viewer.
//
//  Note:
//
//  Programmer: Eric Brugger
//  Creation:   August 4, 2000
//
//  Modifications:
//    Brad Whitlock, Fri Aug 25 11:09:32 PDT 2000
//    I replaced LocalProcess with RemoteProcess since it can now launch
//    local processes.
//
//    Brad Whitlock, Thu Aug 31 14:48:00 PST 2000
//    I made the state objects (PcAttributes) an instance instead of
//    a pointer. This is so observers can be set up before the proxy's
//    remote process is created.
//
//    Eric Brugger, Tue Sep  5 10:21:18 PDT 2000
//    I changed the PlotType enumerated type so that PC_PLOT would be zero.
//    This is currently necessary since it is the only type implemented.
//    Once they are all implemented this will no longer be the case.
//
//    Eric Brugger, Fri Sep 15 11:22:49 PDT 2000
//    I added the methods GetAsliceAttributes and GetGlobalAttributes.
//
//    Eric Brugger, Mon Sep 18 11:42:45 PDT 2000
//    I added MAX_PLOT and MAX_OPERATOR which give a count of the number
//    of plots and operators.
//
//    Brad Whitlock, Tue Sep 19 18:50:18 PST 2000
//    I added a HostProfileList state object.
//
//    Brad Whitlock, Mon Sep 25 11:59:09 PDT 2000
//    I removed PlotType and OperType and put them in include files.
//
//    Brad Whitlock, Thu Sep 28 11:49:47 PDT 2000
//    I added a WriteConfigFile rpc.
//
//    Kathleen Bonnell, Wed Oct 11 08:38:57 PDT 2000
//    I added OnionPeelAttributes.
//
//    Eric Brugger, Wed Oct 25 14:42:06 PDT 2000
//    I removed the argument "name" from the Create method.
//
//    Brad Whitlock, Fri Nov 10 15:01:31 PST 2000
//    I added MaterialAttributes.
//
//    Brad Whitlock, Tue Nov 21 11:25:14 PDT 2000
//    I added a ConnectToMetaDataServer rpc.
//
//    Brad Whitlock, Wed Dec 13 11:06:46 PDT 2000
//    I added MaterialSelectAttributes.
//
//    Hank Childs, Wed Jan 10 11:55:30 PST 2001
//    Add volume attributes.
//
//    Hank Childs, Sun Feb 11 19:00:06 PST 2001
//    Add save window rpc.
//
//    Brad Whitlock, Fri Feb 9 14:20:11 PST 2001
//    Added SaveImageAttributes.
//
//    Brad Whitlock, Fri Feb 16 13:34:34 PST 2001
//    Added ContourAttributes.
//
//    Kathleen Bonnell, Tue Mar  6 10:25:25 PST 2001 
//    Added SurfaceAttributes.
//
//    Eric Brugger, Thu Mar  8 13:00:36 PST 2001
//    Modified to treat plots generically.
//
//    Brad Whitlock, Thu Apr 19 10:54:33 PDT 2001
//    Added methods to handle window iconification.
//
//    Brad Whitlock, Mon Apr 23 09:36:25 PDT 2001
//    Added state objects that can be observed to get error messages and
//    status information.
//
//    Brad Whitlock, Tue Apr 24 10:33:23 PDT 2001
//    Added consecutiveReadZeroes member.
//
//    Brad Whitlock, Mon Apr 30 12:17:31 PDT 2001
//    Added EngineList and StatusAttributes and methods for engine
//    interruption and termination.
//
//    Brad Whitlock, Mon Jun 11 14:18:06 PST 2001
//    Added colortable state object.
//
//    Brad Whitlock, Sun Jun 17 20:07:40 PST 2001
//    Added the AnnotationAttributes object.
//
//    Brad Whitlock, Thu Jun 21 13:00:35 PST 2001
//    Added methods to transfer SIL restrictions.
//
//    Hank Childs, Mon Jul 23 13:43:41 PDT 2001
//    Removed material select.
//
//    Jeremy Meredith, Thu Jul 26 03:11:41 PDT 2001
//    Removed all references to OperType.
//    Added support for real operator plugins.
//
//    Brad Whitlock, Thu Jul 19 15:56:39 PST 2001
//    Added methods to set the view.
//
//    Brad Whitlock, Tue Aug 14 14:52:39 PST 2001
//    Added methods to reset the plot and operator attributes.
//
//    Brad Whitlock, Thu Aug 30 09:52:12 PDT 2001
//    Added methods to set the default annotation attributes.
//
//    Jeremy Meredith, Wed Sep  5 14:06:21 PDT 2001
//    Added plugin manager attributes.
//
//    Brad Whitlock, Tue Sep 4 22:30:11 PST 2001
//    Added appearance attributes.
//
//    Brad Whitlock, Mon Sep 24 11:29:10 PDT 2001
//    Added a method to query the name of the local machine.
//
//    Jeremy Meredith, Fri Sep 28 13:41:57 PDT 2001
//    Added the LoadPlugins method.
//
//    Brad Whitlock, Fri Sep 14 13:45:59 PST 2001
//    Added a light list and RPC's to use the lights.
//
//    Eric Brugger, Mon Nov 19 13:29:49 PST 2001
//    Added animation attributes.
//
//    Brad Whitlock, Wed Sep 19 14:37:18 PST 2001
//    Added RPC's for disabling updates and redrawing the window.
//
//    Kathleen Bonnell, Wed Dec  5 13:42:07 PST 2001
//    Added pick attributes.
//
//    Brad Whitlock, Tue Jan 29 16:16:17 PST 2002
//    Added an RPC to set the window area for the vis windows.
//
//    Brad Whitlock, Wed Feb 20 13:58:29 PST 2002
//    Added printer attributes and RPCs for printing. Also added a method
//    to return the user name.
//
//    Sean Ahern, Tue Apr 16 12:31:39 PDT 2002
//    Added the ability to show/hide all windows.
//
//    Brad Whitlock, Mon May 6 16:32:43 PST 2002
//    Added a bunch of new methods to expose the functions in the
//    viswindow's popup menu.
//
//    Jeremy Meredith, Wed May  8 12:27:29 PDT 2002
//    Added keyframe attributes.
//
//    Hank Childs, Thu May 23 18:36:38 PDT 2002
//    Renamed saveImageAtts to saveWindowAtts.
//
//    Brad Whitlock, Thu Jun 27 16:13:55 PST 2002
//    Added copy window methods.
//
//    Brad Whitlock, Mon Jul 29 15:15:34 PST 2002
//    Added ReOpenDatabase method and ClearCache method.
//
//    Brad Whitlock, Mon Sep 16 12:37:42 PDT 2002
//    I added methods to clear reflines and set view extents. I also
//    added two new state objects.
//
//    Brad Whitlock, Fri Sep 6 13:52:50 PST 2002
//    I added query methods.
//
//    Brad Whitlock, Tue Oct 15 16:24:36 PST 2002
//    I added CloneWindow and CopyPlotsToWindow.
//
//    Jeremy Meredith, Thu Oct 24 16:03:25 PDT 2002
//    Added material options.
//
//    Brad Whitlock, Mon Nov 11 11:45:13 PDT 2002
//    I added methods to lock time and tools.
//
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002 
//    Removed SetPickAttributes, SetDefaultPickAttributes, ResetPickAttributes.
//
//    Eric Brugger, Mon Nov 18 11:58:56 PST 2002
//    Added SetPlotFrameRange and DeletePlotKeyframe.
//
//    Brad Whitlock, Wed Nov 20 14:52:47 PST 2002
//    I changed some of the color table methods.
//
//    Hank Childs, Mon Dec  2 14:17:18 PST 2002
//    Used reference counted pointer with SIL Restriction to be consistent with
//    the rest of the code.
//
//    Eric Brugger, Mon Dec 30 11:15:26 PST 2002
//    Added SetPlotDatabaseState and DeletePlotDatabaseKeyframe.
//
//    Eric Brugger, Fri Jan  3 15:16:07 PST 2003
//    Added ClearViewKeyframes, DeleteViewKeyframe, SetViewKeyframe and
//    ToggleCameraViewMode.
//
//    Brad Whitlock, Thu Dec 19 11:41:48 PDT 2002
//    I added a security key argument to ConnectToMetaDataServer and made
//    changes to interface because ViewerRPC is now automatically generated.
//
//    Brad Whitlock, Mon Jan 13 08:42:12 PDT 2003
//    I added a method to open an mdserver.
//
//    Eric Brugger, Tue Jan 28 12:33:25 PST 2003
//    I added MovePlotKeyframe, MovePlotDatabaseKeyframe and MoveViewKeyframe.
//    I added a plotId argument to SetPlotDatabaseState,
//    DeletePlotDatabaseKeyframe, DeletePlotKeyframe and SetPlotFrameRange.
//
//    Brad Whitlock, Thu Feb 27 11:35:50 PDT 2003
//    I added the animationStopOpcode member.
//
//    Kathleen Bonnell, Tue Mar  4 13:27:11 PST 2003   
//    Added GlobalLineoutAttributes and Set/Get methods. 
//
//    Brad Whitlock, Wed Mar 12 10:45:35 PDT 2003
//    I added the iconifyOpcode member.
//
//    Brad Whitlock, Thu Apr 10 09:27:16 PDT 2003
//    I added PromoteOperator, DemoteOperator, RemoveOperator methods. I also
//    added an overloaded version of SetActivePlots.
//
//    Eric Brugger, Fri Apr 18 15:37:54 PDT 2003
//    I added ToggleMaintainViewMode and deleted ToggleAutoCenterMode.
//
//    Kathleen Bonnell, Thu May 15 10:00:02 PDT 2003  
//    I added ToggleFullFrameMode.
//
//    Brad Whitlock, Thu May 15 13:03:04 PST 2003
//    I added a default timeState argument to the OpenDatabase method.
//
//    Kathleen Bonnell, Tue Jul  1 09:34:37 PDT 2003  
//    Added SetPickAttributes.
//
//    Brad Whitlock, Tue Jul 1 16:45:31 PST 2003
//    I added a method to export color tables.
//
//    Brad Whitlock, Wed Jul 9 11:54:15 PDT 2003
//    I added methods to export and import the viewer's entire state.
//
//    Kathleen Bonnell, Wed Jul 23 17:04:10 PDT 2003
//    Added 'samples' arg to LineQuery and Lineout methods.
//    Added int args to DatabaseQuery.
//    Added overloaded Pick and NodePick methods (accepting doubles).
//
//    Brad Whitlock, Wed Jul 30 14:44:34 PST 2003
//    Added an extra argument to ImportEntireState.
//
//    Eric Brugger, Wed Aug 20 10:44:10 PDT 2003
//    I added GetViewCurveAttributes and SetViewCurve.  I split the view
//    attributes into 2d and 3d parts.
//
//    Brad Whitlock, Fri Aug 29 11:20:32 PDT 2003
//    Added HideToolbars and ShowToolbars.
//
//    Kathleen Bonnell, Thu Sep 11 11:35:08 PDT 2003 
//    Added optional bool arg to AddOperator.
//
//    Brad Whitlock, Wed Oct 15 15:37:20 PST 2003
//    Added optional timeState to ReplaceDatabase.
//
//    Brad Whitlock, Wed Oct 22 12:21:40 PDT 2003
//    Added optional argument to OpenDatabase.
//
//    Kathleen Bonnell, Wed Nov 26 14:17:55 PST 2003 
//    Added ResetPickAttributes. Added optional int args to PointQuery.
//
//    Brad Whitlock, Wed Oct 29 10:31:58 PDT 2003
//    Added new methods to deal with advanced annotations.
//
//    Brad Whitlock, Mon Dec 29 09:17:30 PDT 2003
//    Added methods to set and update the center of rotation.
//
//    Brad Whitlock, Fri Jan 23 09:07:48 PDT 2004
//    I added a list of DatabaseCorrelations and the list of sources. I also
//    renamed a few of the animation methods so they make more sense for
//    multiple time sliders. I added the ActivateDatabase and CheckForNewStates
//    methods. I also added DefineDatabaseCorrelation, AlterDatabaseCorrelation,
//    DeleteDatabaseCorrelation, and CloseDatabase.
//
//    Brad Whitlock, Thu Feb 26 13:37:02 PST 2004
//    Added ClearCacheForAllEngines.
//
//    Jeremy Meredith, Fri Mar 26 10:22:36 PST 2004
//    Added support for simulations.
//
//    Eric Brugger, Mon Mar 29 13:39:53 PST 2004
//    I added ToggleMaintainDataMode.
//
//    Kathleen Bonnell, Wed Mar 31 10:56:30 PST 2004 
//    Added queryOverTimeAtts, bool arg to DatabaseQuery, PointQuery.
//
//    Kathleen Bonnell, Thu Aug  5 08:34:15 PDT 2004 
//    Added ResetLineoutColor.
//
//    Kathleen Bonnell, Wed Aug 18 09:28:51 PDT 2004 
//    Added interactorAtts.
//
//    Jeremy Meredith, Thu Apr 22 14:02:52 PDT 2004
//    Added metaData and GetDatabaseMetaData.  Added silAtts and GetSILAtts.
//
//    Kathleen Bonnell, Wed Dec 15 17:12:47 PST 2004 
//    Added bool arg to DatabaseQuery and PointQuery. 
//
//    Mark C. Miller, Tue Mar  8 17:59:40 PST 2005
//    Added QueryProcessAttributes and GetProcessAttributes
//
//    Jeremy Meredith, Mon Mar 21 08:50:53 PST 2005
//    Added SendSimulationCommand methods.
//
//    Brad Whitlock, Fri Apr 15 11:02:14 PDT 2005
//    Added postponedAction;
//
//    Hank Childs, Wed May 25 10:38:37 PDT 2005
//    Added DBPluginInfo.
//
//    Mark C. Miller, Tue May 31 20:12:42 PDT 2005
//    Added SetTryHarderCyclesTimes
//
//    Brad Whitlock, Tue May 3 16:01:45 PST 2005
//    Added methods to tell the viewer to launch a VisIt client and added
//    objects necessary to allow the viewer to launch the client. I also added
//    clientMethod, a new state object that lets a client execute "methods" on
//    another client. I added a Detach method that lets the client detach from
//    the viewer and leave the other clients running.
//
//    Kathleen Bonnell, Wed Jul 27 15:52:59 PDT 2005
//    Added SuppressQueryOutput.
//
//    Mark C. Miller, Wed Nov 16 10:46:36 PST 2005
//    Added mesh management attributes methods 
//
//    Brad Whitlock, Thu Nov 17 16:37:06 PST 2005
//    Added methods to move and resize windows.
//
//    Hank Childs, Mon Feb 13 21:39:02 PST 2006
//    Added GetConstructDDFAttributes, ConstructDDF.
//
//    Brad Whitlock, Tue Mar 7 16:36:16 PST 2006
//    Added RedoView.
//
//    Kathleen Bonnell, Tue Jun 20 16:02:38 PDT 2006
//    Added UpdatePlotInfoAtts, GetPlotInfoAtts.
//
//    Hank Childs, Mon Jul 10 17:36:41 PDT 2006
//    Add more args to DatabaseQuery.
//
//    Jeremy Meredith, Mon Aug 28 16:55:01 EDT 2006
//    Added ability to force using a specific plugin when opening a file.
//
// ****************************************************************************

class VIEWER_PROXY_API ViewerProxy : public SimpleObserver
{
  public:
    ViewerProxy();
    virtual ~ViewerProxy();

    Connection *GetReadConnection() const;
    Connection *GetWriteConnection() const;
    const std::string &GetLocalHostName() const;
    const std::string &GetLocalUserName() const;
    void ProcessInput();

    void AddArgument(const std::string &arg);
    void Create(int *argc = 0, char ***argv = 0);
    void Close();
    void Detach();
    void LoadPlugins();

    void AddWindow();
    void CloneWindow();
    void DeleteWindow();
    void SetWindowLayout(int layout);
    void SetActiveWindow(int windowId);
    void IconifyAllWindows();
    void DeIconifyAllWindows();
    void ShowAllWindows();
    void HideAllWindows();
    void ClearWindow();
    void ClearAllWindows();
    void SaveWindow();
    void PrintWindow();
    void DisableRedraw();
    void RedrawWindow();
    void ResizeWindow(int win, int w, int h);
    void MoveWindow(int win, int x, int y);
    void MoveAndResizeWindow(int win, int x, int y, int w, int h);
    void HideToolbars(bool forAllWindows = false);
    void ShowToolbars(bool forAllWindows = false);

    void ConnectToMetaDataServer(const std::string &hostName, const stringVector &argv);
    void OpenMDServer(const std::string &hostName, const stringVector &argv);

    void OpenDatabase(const std::string &database, int timeState = 0,
                      bool addDefaultPlots = true,
                      const std::string &forcedFileType = "");
    void CloseDatabase(const std::string &database);
    void ActivateDatabase(const std::string &database);
    void CheckForNewStates(const std::string &database);
    void ReOpenDatabase(const std::string &database, bool forceClose = true);
    void ReplaceDatabase(const std::string &database, int timeState = 0);
    void OverlayDatabase(const std::string &database);
    void ClearCache(const std::string &hostName, const std::string &simName);
    void ClearCacheForAllEngines();
    void UpdateDBPluginInfo(const std::string &hostName);
    void ExportDatabase(void);
    void ConstructDDF(void);

    void CreateDatabaseCorrelation(const std::string &name,
                                   const stringVector &dbs, int method,
                                   int nStates = -1);
    void AlterDatabaseCorrelation(const std::string &name,
                                  const stringVector &dbs, int method,
                                  int nStates = -1);
    void DeleteDatabaseCorrelation(const std::string &name);

    void OpenComputeEngine(const std::string &hostName, const stringVector &argv);
    void CloseComputeEngine(const std::string &hostName, const std::string &simName);
    void InterruptComputeEngine(const std::string &hostName, const std::string &simName);

    void AnimationSetNFrames(int nFrames);
    void AnimationPlay();
    void AnimationReversePlay();
    void AnimationStop();
    void TimeSliderNextState();
    void TimeSliderPreviousState();
    void SetTimeSliderState(int state);
    void SetActiveTimeSlider(const std::string &ts);

    void AddPlot(int type, const std::string &var);
    void SetPlotFrameRange(int plotId, int frame0, int frame1);
    void DeletePlotKeyframe(int plotId, int frame);
    void MovePlotKeyframe(int plotId, int oldFrame, int newFrame);
    void SetPlotDatabaseState(int plotId, int frame, int state);
    void DeletePlotDatabaseKeyframe(int plotId, int frame);
    void MovePlotDatabaseKeyframe(int plotId, int oldFrame, int newFrame);
    void DeleteActivePlots();
    void HideActivePlots();
    void DrawPlots();
    void SetActivePlots(const intVector &activePlotIds);
    void SetActivePlots(const intVector &activePlotIds,
                        const intVector &activeOperatorIds,
                        const intVector &expandedPlots);
    void ChangeActivePlotsVar(const std::string &var);

    void AddOperator(int oper, const bool fromDefault = true);
    void PromoteOperator(int operatorId);
    void DemoteOperator(int operatorId);
    void RemoveOperator(int operatorId);
    void RemoveLastOperator();
    void RemoveAllOperators();

    void SetDefaultPlotOptions(int type);
    void SetPlotOptions(int type);
    void ResetPlotOptions(int type);
    void SetDefaultOperatorOptions(int oper);
    void SetOperatorOptions(int oper);
    void ResetOperatorOptions(int type);

    void SetActiveContinuousColorTable(const std::string &colorTableName);
    void SetActiveDiscreteColorTable(const std::string &colorTableName);
    void DeleteColorTable(const std::string &colorTableName);
    void UpdateColorTable(const std::string &colorTableName);
    void ExportColorTable(const std::string &colorTableName);
    void InvertBackgroundColor();

    void SetCenterOfRotation(double, double, double);
    void ChooseCenterOfRotation();
    void ChooseCenterOfRotation(double, double);
    void SetViewCurve();
    void SetView2D();
    void SetView3D();
    void ClearViewKeyframes();
    void DeleteViewKeyframe(int frame);
    void MoveViewKeyframe(int oldFrame, int newFrame);
    void SetViewKeyframe();
    void ResetView();
    void RecenterView();
    void SetViewExtentsType(int t);
    void ToggleMaintainViewMode();
    void ToggleMaintainDataMode();
    void UndoView();
    void RedoView();
    void ToggleLockViewMode();
    void ToggleLockTime();
    void ToggleLockTools();
    void ToggleSpinMode();
    void ToggleCameraViewMode();
    void ToggleFullFrameMode();

    void SetWindowMode(int mode);
    void ToggleBoundingBoxMode();
    void EnableTool(int tool, bool enabled);

    void CopyViewToWindow(int from, int to);
    void CopyLightingToWindow(int from, int to);
    void CopyAnnotationsToWindow(int from, int to);
    void CopyPlotsToWindow(int from, int to);

    void SetAnnotationAttributes();
    void SetDefaultAnnotationAttributes();
    void ResetAnnotationAttributes();
    void AddAnnotationObject(int annotType);
    void HideActiveAnnotationObjects();
    void DeleteActiveAnnotationObjects();
    void RaiseActiveAnnotationObjects();
    void LowerActiveAnnotationObjects();
    void SetAnnotationObjectOptions();
    void SetDefaultAnnotationObjectList();
    void ResetAnnotationObjectList();

    void SetInteractorAttributes();
    void SetDefaultInteractorAttributes();
    void ResetInteractorAttributes();

    void SetKeyframeAttributes();

    void SetMaterialAttributes();
    void SetDefaultMaterialAttributes();
    void ResetMaterialAttributes();

    void SetLightList();
    void SetDefaultLightList();
    void ResetLightList();

    void SetAnimationAttributes();

    void SetAppearanceAttributes();
    void ProcessExpressions();

    void ClearPickPoints();
    void ClearReferenceLines();

    void SetRenderingAttributes();
    void SetWindowArea(int x, int y, int w, int h);

    void SetGlobalLineoutAttributes();
    void SetPickAttributes();
    void SetDefaultPickAttributes();
    void ResetPickAttributes();
    void ResetPickLetter();

    void ResetLineoutColor();

    void SetQueryOverTimeAttributes();
    void SetDefaultQueryOverTimeAttributes();
    void ResetQueryOverTimeAttributes();

    void SetTryHarderCyclesTimes(int flag);

    void SetMeshManagementAttributes();
    void SetDefaultMeshManagementAttributes();
    void ResetMeshManagementAttributes();

    void WriteConfigFile();
    void ExportEntireState(const std::string &filename);
    void ImportEntireState(const std::string &filename, bool inVisItDir);

    // Methods for dealing with plot SIL restrictions.
    avtSILRestriction_p GetPlotSILRestriction() 
                     { return internalSILRestriction; };
    avtSILRestriction_p GetPlotSILRestriction() const
                     { return new avtSILRestriction(internalSILRestriction); };
    void SetPlotSILRestriction(avtSILRestriction_p newRestriction);
    void SetPlotSILRestriction();

    // Methods for querying
    void SuppressQueryOutput(bool onOff);
    void DatabaseQuery(const std::string &queryName, const stringVector &vars,
                       const bool = false, const int arg1 = 0, const int arg2 = 0,
                       const bool = false, const double darg1 = 0., const double darg2 = 0.);
    void PointQuery(const std::string &queryName, const double pt[3],
                    const stringVector &vars, const bool = false,
                    const int arg1 = -1, const int arg2 = -1,
                    const bool = false);
    void LineQuery(const std::string &queryName, const double pt1[3],
                   const double pt2[3], const stringVector &vars,
                   const int samples);
    void Pick(int x, int y, const stringVector &vars);
    void Pick(double xyz[3], const stringVector &vars);
    void NodePick(int x, int y, const stringVector &vars);
    void NodePick(double xyz[3], const stringVector &vars);
    void Lineout(const double p0[3], const double p1[3],
                 const stringVector &vars, const int samples);

    void QueryProcessAttributes(const std::string componentName,
                                const std::string engineHostName,
                                const std::string engineDbName);

    void SendSimulationCommand(const std::string &hostName,
                               const std::string &simName,
                               const std::string &command);

    void SendSimulationCommand(const std::string &hostName,
                               const std::string &simName,
                               const std::string &command,
                               const std::string &argument);

    void OpenClient(const std::string &clientName, 
                    const std::string &program,
                    const stringVector &args);
    int MethodRequestHasRequiredInformation() const;

    // Methods for returning pointers to state obects.
    AnimationAttributes        *GetAnimationAttributes() const 
                                    {return animationAtts;};
    AnnotationAttributes       *GetAnnotationAttributes() const 
                                    {return annotationAtts;};
    AppearanceAttributes       *GetAppearanceAttributes() const 
                                    {return appearanceAtts;};
    ColorTableAttributes       *GetColorTableAttributes() const 
                                    {return colorTableAtts;};
    ConstructDDFAttributes     *GetConstructDDFAttributes() const
                                    {return constructDDFAtts;}
    DatabaseCorrelationList    *GetDatabaseCorrelationList() const
                                    {return correlationList; };
    DBPluginInfoAttributes     *GetDBPluginInfoAttributes() const
                                    {return dbPluginInfoAtts;}
    ExportDBAttributes         *GetExportDBAttributes() const
                                    {return exportDBAtts;}
    EngineList                 *GetEngineList() const 
                                    {return engineList;};
    ExpressionList             *GetExpressionList() const 
                                    {return exprList;};
    GlobalAttributes           *GetGlobalAttributes() const 
                                    {return globalAtts;};
    HostProfileList            *GetHostProfileList() const 
                                    {return hostProfiles;};
    InteractorAttributes       *GetInteractorAttributes() const 
                                    {return interactorAtts;};
    KeyframeAttributes         *GetKeyframeAttributes() const
                                    {return keyframeAtts;}
    LightList                  *GetLightList() const 
                                    {return lightList; };
    MessageAttributes          *GetMessageAttributes() const 
                                    {return messageAtts;};
    AttributeSubject           *GetOperatorAttributes(int type) const;
    PickAttributes             *GetPickAttributes() const 
                                    {return pickAtts;};
    QueryAttributes            *GetQueryAttributes() const 
                                    {return queryAtts;};
    AttributeSubject           *GetPlotAttributes(int type) const;
    PlotList                   *GetPlotList() const 
                                    {return plotList;};
    PluginManagerAttributes    *GetPluginManagerAttributes() const 
                                    {return pluginManagerAttributes;};
    PrinterAttributes          *GetPrinterAttributes() const 
                                    {return printerAtts;};
    SaveWindowAttributes       *GetSaveWindowAttributes() const 
                                    {return saveWindowAtts;};
    SILRestrictionAttributes   *GetSILRestrictionAttributes() const 
                                    {return silRestrictionAtts;};
    StatusAttributes           *GetStatusAttributes() const 
                                    {return statusAtts;};
    SyncAttributes             *GetSyncAttributes() const
                                    {return syncAtts; };
    ViewCurveAttributes        *GetViewCurveAttributes() const 
                                    {return viewCurveAttributes;};
    View2DAttributes           *GetView2DAttributes() const 
                                    {return view2DAttributes;};
    View3DAttributes           *GetView3DAttributes() const 
                                    {return view3DAttributes;};
    WindowInformation          *GetWindowInformation() const
                                    {return windowInfo; };
    RenderingAttributes        *GetRenderingAttributes() const
                                    {return renderAtts; };
    QueryList                  *GetQueryList() const
                                    {return queryList; };
    MaterialAttributes         *GetMaterialAttributes() const
                                    {return materialAtts;}
    GlobalLineoutAttributes    *GetGlobalLineoutAttributes() const 
                                    {return globalLineoutAtts;};
    AnnotationObjectList       *GetAnnotationObjectList() const
                                    {return annotationObjectList; };
    QueryOverTimeAttributes    *GetQueryOverTimeAttributes() const 
                                    {return queryOverTimeAtts;};
    avtDatabaseMetaData        *GetDatabaseMetaData() const
                                    {return metaData; }
    SILAttributes              *GetSILAtts() const
                                    {return silAtts; }
    ProcessAttributes          *GetProcessAttributes() const
                                    {return procAtts; }
    ClientMethod               *GetClientMethod() const
                                    {return clientMethod; }
    ClientInformation          *GetClientInformation() const
                                    {return clientInformation; }
    const ClientInformationList *GetClientInformationList() const
                                    {return clientInformationList; }
    MovieAttributes            *GetMovieAttributes() const
                                    {return movieAtts; }
    MeshManagementAttributes   *GetMeshManagementAttributes() const
                                    {return meshManagementAtts;}
    ViewerRPC                  *GetLogRPC() const
                                    {return logRPC; }
    PlotInfoAttributes         *GetPlotInfoAtts() const
                                    {return plotInfoAtts; }
    void                        UpdatePlotInfoAtts(int winId=-1, int plotId=-1);

    // Don't use this method unless absolutely necessary.
    void SetXferUpdate(bool val);

  protected:
    virtual void Update(Subject *subj);
  private:
    void ConnectXfer();

    RemoteProcess              *viewer;
    ParentProcess              *viewerP;
    Xfer                       *xfer;
    ViewerRPC                  *viewerRPC;

    int                        nPlots;
    int                        nOperators;
    int                        animationStopOpcode;
    int                        iconifyOpcode;

    // State objects
    PostponedAction            *postponedAction;
    SyncAttributes             *syncAtts;
    GlobalAttributes           *globalAtts;
    DatabaseCorrelationList    *correlationList;
    DBPluginInfoAttributes     *dbPluginInfoAtts;
    ExportDBAttributes         *exportDBAtts;
    ConstructDDFAttributes     *constructDDFAtts;
    PlotList                   *plotList;
    ColorTableAttributes       *colorTableAtts;
    ExpressionList             *exprList;
    HostProfileList            *hostProfiles;
    InteractorAttributes       *interactorAtts;
    MessageAttributes          *messageAtts;
    SaveWindowAttributes       *saveWindowAtts;
    StatusAttributes           *statusAtts;
    EngineList                 *engineList;
    AnnotationAttributes       *annotationAtts;
    SILRestrictionAttributes   *silRestrictionAtts;
    ViewCurveAttributes        *viewCurveAttributes;
    View2DAttributes           *view2DAttributes;
    View3DAttributes           *view3DAttributes;
    LightList                  *lightList;
    MaterialAttributes         *materialAtts;
    AnimationAttributes        *animationAtts;
    PluginManagerAttributes    *pluginManagerAttributes;
    AppearanceAttributes       *appearanceAtts;
    PickAttributes             *pickAtts;
    PrinterAttributes          *printerAtts;
    KeyframeAttributes         *keyframeAtts;
    WindowInformation          *windowInfo;
    RenderingAttributes        *renderAtts;
    QueryList                  *queryList;
    QueryAttributes            *queryAtts;
    GlobalLineoutAttributes    *globalLineoutAtts;
    AnnotationObjectList       *annotationObjectList;
    QueryOverTimeAttributes    *queryOverTimeAtts;
    avtDatabaseMetaData        *metaData;
    SILAttributes              *silAtts;
    ProcessAttributes          *procAtts;
    ClientMethod               *clientMethod;
    ClientInformation          *clientInformation;
    ClientInformationList      *clientInformationList;
    MovieAttributes            *movieAtts;
    MeshManagementAttributes   *meshManagementAtts;
    ViewerRPC                  *logRPC;
    PlotInfoAttributes         *plotInfoAtts;

    AttributeSubject           **plotAtts;
    AttributeSubject           **operatorAtts;

    // Extra command line arguments to pass to the viewer.
    stringVector               argv;

    // Used to store the sil restriction in avt format.
    avtSILRestriction_p        internalSILRestriction;
};

#endif
