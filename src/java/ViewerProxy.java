// ****************************************************************************
//
// Copyright (c) 2000 - 2006, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ****************************************************************************

package llnl.visit;

import java.lang.ArrayIndexOutOfBoundsException;
import java.util.Vector;
import java.util.prefs.Preferences;
import java.util.prefs.BackingStoreException;

// ****************************************************************************
// Class: ViewerProxy
//
// Purpose:
//   This is the main class that users of the Java VisIt Interface need to
//   use in order to control VisIt from Java. This class provides methods to
//   launch VisIt's viewer tell it to do things.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Aug 8 12:43:40 PDT 2002
//
// Modifications:
//   Brad Whitlock, Mon Sep 16 12:46:02 PDT 2002
//   I added a method to clear reflines and a method to set the view extents.
//   I also added the windowInfo and renderAtts members.
//
//   Brad Whitlock, Tue Sep 24 07:22:27 PDT 2002
//   I fixed the Create method so it blocks until the plugins are loaded.
//   This prevents an error where an exception is thrown when looking up the
//   name of an exception before the plugins are loaded. Observers can still
//   attach to the subjects in the ViewerProxy but in order to get initial
//   state from the viewer, they should be connected before the call to the
//   Create method.
//
//   Brad Whitlock, Tue Oct 15 16:22:42 PST 2002
//   I added methods to clone windows and copy plots from one window to
//   another window.
//
//   Jeremy Meredith, Thu Oct 24 16:15:11 PDT 2002
//   Added material options.
//
//   Brad Whitlock, Mon Nov 11 11:43:28 PDT 2002
//   I added methods to lock tools and lock time.
//
//   Brad Whitlock, Thu Dec 12 10:22:22 PDT 2002
//   I fixed some color table issues.
//
//   Brad Whitlock, Fri Dec 13 17:03:00 PST 2002
//   I added keyframing routines.
//
//   Brad Whitlock, Fri Dec 27 15:45:38 PST 2002
//   I made it use a new automatically generated version of ViewerRPC so it
//   is trivial to keep the code up to date with the C++ version. I also
//   added new methods that were added while I was making this update.
//
//   Brad Whitlock, Thu Mar 20 10:30:13 PDT 2003
//   I added GlobalLineoutAttributes.
//
//   Brad Whitlock, Thu Apr 10 09:27:16 PDT 2003
//   I added PromoteOperator, DemoteOperator, RemoveOperator methods. I also
//   added an override of SetActivePlots.
//
//   Brad Whitlock, Mon May 5 14:13:33 PST 2003
//   I made it use the regenerated ViewerRPC.java class.
//
//   Brad Whitlock, Thu May 15 13:05:09 PST 2003
//   I added a new OpenDatabase method that lets us open databases at late
//   time states.
//
//   Kathleen Bonnell, Tue Jul  1 10:11:37 PDT 2003
//   I added SetPickAttributes, SetGlobalLineoutAttributes. 
//
//   Brad Whitlock, Tue Jul 1 16:49:48 PST 2003
//   I added a method to tell the viewer to export a colortable.
//
//   Brad Whitlock, Wed Jul 9 12:30:58 PDT 2003
//   I added methods to tell the viewer to export and import its entire state.
//
//   Eric Brugger, Wed Aug 27 09:10:35 PDT 2003
//   I added viewCurveAttributes.  I split the view attributes into 2d
//   and 3d parts.
//
//   Brad Whitlock, Fri Aug 29 11:24:11 PDT 2003
//   I added HideToolbars and ShowToolbars.
//
//   Brad Whitlock, Wed Oct 15 15:39:02 PST 2003
//   I added a timeState argument to the ReplaceDatabase method.
//
//   Brad Whitlock, Wed Oct 22 12:25:43 PDT 2003
//   I added another overloaded OpenDatabase method.
//
//   Brad Whitlock, Wed Oct 29 10:41:21 PDT 2003
//   I added new RPCs to deal with setting up advanced annotations.
//
//   Kathleen Bonnell, Wed Nov 26 14:16:53 PST 2003 
//   Added ResetPickAttributes method. Added overloaded PointQuery and
//   DatabaseQuery methods to handle optional arguments.
//
//   Kathleen Bonnell, Wed Dec 17 15:19:46 PST 2003 
//   Added SetDefaultPickAttributes, ResetPickLetter. 
//
//   Brad Whitlock, Mon Dec 29 09:24:20 PDT 2003
//   Added methods for setting and choosing the center of rotation.
//
//   Brad Whitlock, Thu Jan 8 14:19:43 PST 2004
//   I fixed a typo that prevented it from building.
//
//   Brad Whitlock, Fri Jan 23 14:15:49 PST 2004
//   I added DatabaseCorrelationList and the CloseDatabase rpc.
//
//   Brad Whitlock, Thu Feb 26 13:38:26 PST 2004
//   I added ClearCacheForAllEngines.
//
//   Jeremy Meredith, Tue Mar 30 17:36:24 PST 2004
//   Added support for a simulation name in addition to simply a host name
//   when specifying engines (i.e. ClearCache and CloseComputeEngine).
//
//   Kathleen Bonnell, Wed Mar 31 07:35:25 PST 2004 
//   Added support for QueryOverTimeAttributes.
//
//   Brad Whitlock, Mon Jul 26 15:57:32 PST 2004
//   Added CheckForNewStates.
//
//   Kathleen Bonnell, Wed Aug 18 09:19:25 PDT 2004 
//   Added support for InteractorAttributes.
//
//   Brad Whitlock, Fri Jan 7 17:00:00 PST 2005
//   I hooked up the InteractorAttributes to xfer and added a silAtts
//   member to get the sil attributes for simulations.
//
//   Brad Whitlock, Fri Apr 8 11:57:46 PDT 2005
//   Overloaded AddOperator.
//
//   Brad Whitlock, Fri Apr 15 10:32:35 PDT 2005
//   Added PostponedAction.
//
//   Brad Whitlock, Thu May 12 13:46:37 PST 2005
//   Added ProcessAttributes.
//
//   Brad Whitlock, Wed May 4 09:53:34 PDT 2005
//   Changed the order of some state objects in xfer and added support for
//   reverse launching. I also added OpenClient, AssistOpenClient, and the
//   new clientMethod, clientInformation, clientInformationList, and
//   movieAttributes members.
//
//   Brad Whitlock, Thu Jun 2 11:15:40 PDT 2005
//   Added SetTryHarderCyclesTimes.
//
//   Brad Whitlock, Mon Jun 6 09:39:33 PDT 2005
//   Added some code to yield to other threads so we don't tie up the CPU
//   in some of our blocking loops. I also added code to try and determine
//   VisIt's installation location using java.util.prefs. The Windows
//   installer for VisIt pokes the right values into the registry so Java
//   will know where to look for VisIt. On UNIX, we use the default binary
//   path unless the user overrides it. I also added GetDataPath so it is
//   easier to write test programs that know where to find default data files.
//
//   Brad Whitlock, Thu Jun 16 17:24:58 PST 2005
//   I made the code that waits for the plugins to get loaded use the Yielder
//   class, which makes the thread sleep if yield is not sufficient to make
//   the thread reduce its CPU usage. I also changed things so the Create
//   method never turns off xfer's reading thread. This means that the class
//   will now always read from the viewer unless you call StopProcessing.
//
//   Brad Whitlock, Thu Jul 14 12:25:11 PDT 2005
//   I made SetPlotOptions and SetOperatorOptions complain if they are given
//   invalid indices or names.
//
//   Brad Whitlock, Thu Oct 27 10:45:59 PDT 2005
//   I made it possible to skip plugins that did not load while still allowing
//   plugins that did load to get hooked up to xfer.
//
//   Brad Whitlock, Thu Nov 17 16:36:17 PST 2005
//   I added methods to move and resize windows.
//
//   Brad Whitlock, Tue Nov 29 12:13:41 PDT 2005
//   I added mesh management attributes.
//
//   Brad Whitlock, Tue Jan 10 16:29:28 PST 2006
//   Added another ViewerRPC object for logging.
//
//   Brad Whitlock, Tue Mar 7 16:43:46 PST 2006
//   Added RedoView.
//
//   Brad Whitlock, Thu Mar 16 16:23:50 PST 2006
//   Added DDF attributes.
//
//   Brad Whitlock, Fri Sep 22 14:17:25 PST 2006
//   Added PlotInfoAttributes.
//
// ****************************************************************************

public class ViewerProxy implements SimpleObserver
{
    public ViewerProxy()
    {
        viewer = new RemoteProcess("visit");
        viewer.AddArgument("-viewer");

        // Make viewer RPC's synchronous by default.
        syncCount = 0;
        synchronous = true;
        pluginsLoaded = false;
        doUpdate = true;
        verbose = false;

        xfer = new Xfer();
        messageObserver = new MessageObserver();
        syncNotifier = new SyncNotifier();

        // State objects
        rpc = new ViewerRPC();
        postponedAction = new PostponedAction();
        syncAtts = new SyncAttributes();
        appearanceAtts = new AppearanceAttributes();
        pluginAtts = new PluginManagerAttributes();
        globalAtts = new GlobalAttributes();
        correlationList = new DatabaseCorrelationList();
        plotList = new PlotList();
        hostProfiles = new HostProfileList();
        messageAtts = new MessageAttributes();
        saveAtts = new SaveWindowAttributes();
        statusAtts = new StatusAttributes();
        engineList = new EngineList();
        colorTableAtts = new ColorTableAttributes();
        expressionList = new ExpressionList();
        annotationAtts = new AnnotationAttributes();
        silRestrictionAtts = new SILRestrictionAttributes();
        viewCurve = new ViewCurveAttributes();
        view2D = new View2DAttributes();
        view3D = new View3DAttributes();
        lightList = new LightList();
        materialAtts = new MaterialAttributes();
        animationAtts = new AnimationAttributes();
        pickAtts = new PickAttributes();
        printerAtts = new PrinterAttributes();
        windowInfo = new WindowInformation();
        renderAtts = new RenderingAttributes();
        keyframeAtts = new KeyframeAttributes();
        queryList = new QueryList();
        queryAtts = new QueryAttributes();
        globalLineoutAtts = new GlobalLineoutAttributes();
        annotationObjectList = new AnnotationObjectList();
        queryOverTimeAtts = new QueryOverTimeAttributes();
        interactorAtts = new InteractorAttributes();
        silAtts = new SILRestrictionAttributes();
        processAtts = new ProcessAttributes();
        dbPluginInfoAtts = new DBPluginInfoAttributes();
        exportDBAtts = new ExportDBAttributes();
        constructDDFAtts = new ConstructDDFAttributes();
        clientMethod = new ClientMethod();
        clientInformation = new ClientInformation();
        clientInformationList = new ClientInformationList();
        movieAtts = new MovieAttributes();
        meshManagementAtts = new MeshManagementAttributes();
        logRPC = new ViewerRPC();
        plotInfoAtts = new PlotInfoAttributes();

        // Create the plugin managers.
        plotPlugins = new PluginManager("plot");
        operatorPlugins = new PluginManager("operator");

        // The VISITHOME Java preference will not have been set on
        // UNIX, set the default paths to where we install in /usr/gapps.
        SetBinPath("/usr/gapps/visit/bin");
        dataPath = new String("/usr/gapps/visit/data/");

        // Now that the default paths have been set, try and override them
        // with information from the Java preferences.
        try
        {
            String visitHomeKey = new String("/llnl/visit");           
            if(Preferences.systemRoot().nodeExists(visitHomeKey))
            {
                // Get the node
                Preferences visitNode = Preferences.systemRoot().node(visitHomeKey);
                String defaultValue = new String("");
                String visitPath = visitNode.get("VISITHOME", defaultValue);
                if(!visitPath.equals(defaultValue))
                {
                    PrintMessage("Setting visit path to: " + visitPath);
                    SetBinPath(visitPath);
                    dataPath = new String(visitPath + "\\data\\");
                }
            }
        }
        catch(BackingStoreException e)
        {
            PrintMessage("Cannot access preferences");
        }
    }

    // Close the viewer when the object is garbage collected.
    protected void finalize() throws Throwable
    {
        Close();
        super.finalize();
    }

    // Sets the location of the visit binary.
    public void SetBinPath(String path)
    {
        viewer.SetBinPath(path);
    }

    public String GetDataPath()
    {
        return dataPath;
    }

    // Adds extra arguments to the viewer's command line.
    public void AddArgument(String arg)
    {
        viewer.AddArgument(arg);
    }

    public synchronized void PrintMessage(String msg)
    {
        if(verbose)
            System.out.println(msg);
    }

    // Launches the viewer and sets up state objects.
    public boolean Create(int port)
    {
        boolean retval = false;

        if(viewer.Open(port))
        {
            retval = true;

            PrintMessage("Adding state objects.");

            // Set up xfer and the RPC's
            xfer.SetRemoteProcess(viewer);
            xfer.Add(rpc);
            xfer.Add(postponedAction);
            xfer.Add(syncAtts);
            xfer.Add(messageAtts);
            xfer.Add(statusAtts);
            // For simulations. Note that the metadata is being
            // added as a dummy so Xfer can skip its bytes when
            // messages come in because we don't have an object
            // for it yet.
            xfer.AddDummy(); // metadata
            xfer.Add(silAtts);
            xfer.Add(dbPluginInfoAtts);
            xfer.Add(exportDBAtts);
            xfer.Add(constructDDFAtts);
            xfer.Add(clientMethod);
            xfer.Add(clientInformation);
            xfer.Add(clientInformationList);
            xfer.Add(plotInfoAtts);

            xfer.Add(pluginAtts);
            xfer.Add(appearanceAtts);
            xfer.Add(globalAtts);
            xfer.Add(correlationList);
            xfer.Add(plotList);
            xfer.Add(hostProfiles);
            xfer.Add(saveAtts);
            xfer.Add(engineList);
            xfer.Add(colorTableAtts);
            xfer.Add(expressionList);
            xfer.Add(annotationAtts);
            xfer.Add(silRestrictionAtts);
            xfer.Add(viewCurve);
            xfer.Add(view2D);
            xfer.Add(view3D);
            xfer.Add(lightList);
            xfer.Add(animationAtts);
            xfer.Add(pickAtts);
            xfer.Add(printerAtts);
            xfer.Add(windowInfo);
            xfer.Add(renderAtts);
            xfer.Add(keyframeAtts);
            xfer.Add(queryList);
            xfer.Add(queryAtts);
	    xfer.Add(materialAtts);
            xfer.Add(globalLineoutAtts);
            xfer.Add(annotationObjectList);
            xfer.Add(queryOverTimeAtts);
            xfer.Add(interactorAtts);
            xfer.Add(processAtts);
            xfer.Add(movieAtts);
            xfer.Add(meshManagementAtts);
            xfer.Add(logRPC);

            // hook up the message observer.
            messageObserver.Attach(messageAtts);

            // Hook up the syncAtts notifier.
            syncAtts.Attach(syncNotifier);

            // Hook up this object to the plugin atts.
            pluginAtts.Attach(this);

            // Start reading input from the viewer so we can load the
            // plugins that the viewer has loaded.
            StartProcessing();
            Yielder y = new Yielder(1000);
            boolean errFlag = false;
            while(!errFlag && !pluginsLoaded)
            {
                // Yield to other threads so we don't swamp the CPU with
                // zero work.
                try
                {
                    y.yield();
                }
                catch(java.lang.InterruptedException e)
                {
                    errFlag = true;
                }
            }
        }

        return retval;
    }

    // Starts automatic processing of information from the viewer.
    public void StartProcessing()
    {
        PrintMessage("Starting to read from input from viewer.");
        xfer.StartProcessing();
    }

    // Stops automatic processing of information from the viewer.
    public void StopProcessing()
    {
        PrintMessage("Stopped reading input from viewer.");
        xfer.StopProcessing();
    }

    // Sets the synchronous flag which determines if there is a synchronization
    // step after sending a viewer rpc.
    public void SetSynchronous(boolean val)
    {
        synchronous = val;
    }

    // Sets a flag that determines if messages from the viewer are printed to
    // the console.
    public void SetVerbose(boolean val)
    {
        verbose = val;
        messageObserver.SetVerbose(val);
        plotPlugins.SetVerbose(val);
        operatorPlugins.SetVerbose(val);
        viewer.SetVerbose(val);
        xfer.SetVerbose(val);
        syncNotifier.SetVerbose(val);
    }

//
// Viewer RPC's
//
    public boolean Close()
    {
        PrintMessage("Closing the viewer.");
        xfer.StopProcessing();
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLOSERPC);
        rpc.Notify();
        viewer.Close();
        return true;
    }

    public boolean AddWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CloneWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLONEWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeleteWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DELETEWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetWindowLayout(int layout)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETWINDOWLAYOUTRPC);
        rpc.SetWindowLayout(layout);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActiveWindow(int windowId)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETACTIVEWINDOWRPC);
        rpc.SetWindowId(windowId);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean IconifyAllWindows()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ICONIFYALLWINDOWSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeIconifyAllWindows()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DEICONIFYALLWINDOWSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public void SetWindowArea(int x, int y, int w, int h)
    {
        String tmp = new String(x + "x" + y + "+" + x + "+" + y);
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETWINDOWAREARPC);
        rpc.SetWindowArea(tmp);
        rpc.Notify();
    }

    public void ShowAllWindows()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SHOWALLWINDOWSRPC);
        rpc.Notify();
    }

    public boolean HideAllWindows()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_HIDEALLWINDOWSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearAllWindows()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARALLWINDOWSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SaveWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SAVEWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean PrintWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_PRINTWINDOWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DisableRedraw()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DISABLEREDRAWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RedrawWindow()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REDRAWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OpenDatabase(String database)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OPENDATABASERPC);
        rpc.SetDatabase(database);
        rpc.SetIntArg1(0);
        rpc.SetBoolFlag(true);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OpenDatabase(String database, int timeState)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OPENDATABASERPC);
        rpc.SetDatabase(database);
        rpc.SetIntArg1(timeState);
        rpc.SetBoolFlag(true);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OpenDatabase(String database, int timeState,
                                boolean addDefaultPlots)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OPENDATABASERPC);
        rpc.SetDatabase(database);
        rpc.SetIntArg1(timeState);
        rpc.SetBoolFlag(addDefaultPlots);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CloseDatabase(String database)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLOSEDATABASERPC);
        rpc.SetDatabase(database);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ReOpenDatabase(String database)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REOPENDATABASERPC);
        rpc.SetDatabase(database);
        rpc.SetWindowLayout(1);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CheckForNewStates(String database)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CHECKFORNEWSTATESRPC);
        rpc.SetDatabase(database);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ReplaceDatabase(String database, int timeState)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REPLACEDATABASERPC);
        rpc.SetDatabase(database);
        rpc.SetIntArg1(timeState);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OverlayDatabase(String database)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OVERLAYDATABASERPC);
        rpc.SetDatabase(database);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearCache(String hostName, String simName)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARCACHERPC);
        rpc.SetProgramHost(hostName);
        rpc.SetProgramSim(simName);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearCacheForAllEngines()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARCACHEFORALLENGINESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OpenComputeEngine(String hostName, Vector argv)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OPENCOMPUTEENGINERPC);
        rpc.SetProgramHost(hostName);
        rpc.SetProgramOptions(argv);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CloseComputeEngine(String hostName, String simName)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLOSECOMPUTEENGINERPC);
        rpc.SetProgramHost(hostName);
        rpc.SetProgramSim(simName);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public void InterruptComputeEngine()
    {
        xfer.SendInterruption();
    }

    public boolean AnimationSetNFrames(int nFrames)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ANIMATIONSETNFRAMESRPC);
        rpc.SetNFrames(nFrames);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AnimationPlay()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ANIMATIONPLAYRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AnimationReversePlay()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ANIMATIONREVERSEPLAYRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AnimationStop()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ANIMATIONSTOPRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean TimeSliderNextState()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TIMESLIDERNEXTSTATERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean TimeSliderPreviousState()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TIMESLIDERPREVIOUSSTATERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetTimeSliderState(int state)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETTIMESLIDERSTATERPC);
        rpc.SetStateNumber(state);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AddPlot(int type, String var)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDPLOTRPC);
        rpc.SetPlotType(type);
        rpc.SetVariable(var);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AddPlot(String plotName, String var)
    {
        boolean retval = false;
        int type = plotPlugins.IndexFromName(plotName);
        if(type > -1)
        {
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDPLOTRPC);
            rpc.SetPlotType(type);
            rpc.SetVariable(var);
            rpc.Notify();
            retval = synchronous ? Synchronize() : true;
        }

        return retval;
    }

    public boolean SetPlotFrameRange(int frame0, int frame1)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETPLOTFRAMERANGERPC);
        rpc.SetFrameRange(frame0, frame1);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeletePlotKeyframe(int frame)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DELETEPLOTKEYFRAMERPC);
        rpc.SetFrame(frame);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeleteActivePlots()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DELETEACTIVEPLOTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean HideActivePlots()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_HIDEACTIVEPLOTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DrawPlots()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DRAWPLOTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActivePlots(Vector activePlots, Vector activeOperators,
                                  Vector expandedPlots)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETACTIVEPLOTSRPC);
        rpc.SetActivePlotIds(activePlots);
        rpc.SetActiveOperatorIds(activeOperators);
        rpc.SetExpandedPlotIds(expandedPlots);
        rpc.SetBoolFlag(true);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActivePlots(Vector ids)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETACTIVEPLOTSRPC);
        rpc.SetActivePlotIds(ids);
        rpc.SetBoolFlag(false);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActivePlots(int[] ids)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETACTIVEPLOTSRPC);
        Vector iv = new Vector();
        for(int i = 0; i < ids.length; ++i)
            iv.addElement(new Integer(ids[i]));
        rpc.SetActivePlotIds(iv);
        rpc.SetBoolFlag(false);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActivePlot(int index)
    {
        int[] ids = new int[1];
        ids[0] = index;
        return SetActivePlots(ids);
    }

    public boolean ChangeActivePlotsVar(String var)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CHANGEACTIVEPLOTSVARRPC);
        rpc.SetVariable(var);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AddOperator(int oper)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDOPERATORRPC);
        rpc.SetOperatorType(oper);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AddOperator(String type)
    {
        int oper = operatorPlugins.IndexFromName(type);
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDOPERATORRPC);
        rpc.SetOperatorType(oper);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean PromoteOperator(int operatorId)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_PROMOTEOPERATORRPC);
        rpc.SetOperatorType(operatorId);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DemoteOperator(int operatorId)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DEMOTEOPERATORRPC);
        rpc.SetOperatorType(operatorId);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RemoveOperator(int operatorId)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REMOVEOPERATORRPC);
        rpc.SetOperatorType(operatorId);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RemoveLastOperator()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REMOVELASTOPERATORRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RemoveAllOperators()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REMOVEALLOPERATORSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultPlotOptions(int type)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTPLOTOPTIONSRPC);
        rpc.SetPlotType(type);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetPlotOptions(int type)
    {
        boolean retval = false;
        if(type >= 0)
        {
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETPLOTOPTIONSRPC);
            rpc.SetPlotType(type);
            rpc.Notify();
            retval = synchronous ? Synchronize() : true;
        }
        else
        {
            PrintMessage("SetPlotOptions: " + type + 
                         " is an invalid plot index.");
        }

        return retval;
    }

    public boolean SetPlotOptions(String type)
    {
        return SetPlotOptions(plotPlugins.IndexFromName(type));
    }

    public boolean ResetPlotOptions(int type)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETPLOTOPTIONSRPC);
        rpc.SetPlotType(type);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultOperatorOptions(int type)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTOPERATOROPTIONSRPC);
        rpc.SetOperatorType(type);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetOperatorOptions(int type)
    {
        boolean retval = false;
        if(type >= 0)
        {
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETOPERATOROPTIONSRPC);
            rpc.SetOperatorType(type);
            rpc.Notify();
            retval = synchronous ? Synchronize() : true;
        }
        else
        {
            PrintMessage("SetPlotOptions: " + type + 
                         " is an invalid operator index.");
        }

        return retval;
    }

    public boolean SetOperatorOptions(String type)
    {
        return SetOperatorOptions(operatorPlugins.IndexFromName(type));
    }

    public boolean ResetOperatorOptions(int type)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETOPERATOROPTIONSRPC);
        rpc.SetOperatorType(type);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetActiveContinuousColorTable(String colorTableName)
    {
        if(colorTableAtts.GetColorTableIndex(colorTableName) != -1)
        {
            colorTableAtts.SetActiveContinuous(colorTableName);
            colorTableAtts.Notify();

            // Update the color table. This has the effect of making all plots
            // use the default color table update to use the new active color
            // table.
            boolean sync = synchronous;
            synchronous = false;
            UpdateColorTable(colorTableName);
            synchronous = sync;
        }

        return synchronous ? Synchronize() : true;
    }

    public boolean SetActiveDiscreteColorTable(String colorTableName)
    {
        if(colorTableAtts.GetColorTableIndex(colorTableName) != -1)
        {
            colorTableAtts.SetActiveDiscrete(colorTableName);
            colorTableAtts.Notify();
        }

        return synchronous ? Synchronize() : true;
    }

    public boolean DeleteColorTable(String colorTableName)
    {
        // If it's a valid color table name, make it active.
        int index = colorTableAtts.GetColorTableIndex(colorTableName);
        if(index != -1)
        {
            // Remove the color table from the list and update.
            colorTableAtts.RemoveColorTable(index);
            colorTableAtts.Notify();

            // Update the color table. The specified color table will no
            // longer exist in the list of color tables so all plots that used
            // that color table will have their color tables changed to something
            // else.
            boolean sync = synchronous;
            synchronous = false;
            UpdateColorTable(colorTableName);
            synchronous = sync;
        }

        return synchronous ? Synchronize() : true;
    }

    public boolean UpdateColorTable(String colorTableName)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_UPDATECOLORTABLERPC);
        rpc.SetColorTableName(colorTableName);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ExportColorTable(String colorTableName)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_EXPORTCOLORTABLERPC);
        rpc.SetColorTableName(colorTableName);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean InvertBackgroundColor()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_INVERTBACKGROUNDRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean WriteConfigFile()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_WRITECONFIGFILERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ExportEntireState(String filename)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_EXPORTENTIRESTATERPC);
        rpc.SetVariable(filename);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ImportEntireState(String filename, boolean inVisItDir)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_IMPORTENTIRESTATERPC);
        rpc.SetVariable(filename);
        rpc.SetBoolFlag(inVisItDir);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetRenderingAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETRENDERINGATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }


    public boolean SetCenterOfRotation(double[] pt)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETCENTEROFROTATIONRPC);
        rpc.SetQueryPoint1(pt);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ChooseCenterOfRotation()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CHOOSECENTEROFROTATIONRPC);
        rpc.SetBoolFlag(false);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ChooseCenterOfRotation(double[] pt)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CHOOSECENTEROFROTATIONRPC);
        rpc.SetBoolFlag(true);
        rpc.SetQueryPoint1(pt);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetViewCurve()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETVIEWCURVERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetView2D()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETVIEW2DRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetView3D()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETVIEW3DRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearViewKeyframes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARVIEWKEYFRAMESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeleteViewKeyframe(int frame)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DELETEVIEWKEYFRAMERPC);
        rpc.SetFrame(frame);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetViewKeyframe()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETVIEWKEYFRAMERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetView()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETVIEWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RecenterView()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RECENTERVIEWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleMaintainViewMode()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLEMAINTAINVIEWMODERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleSpinMode()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLESPINMODERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleCameraViewMode()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLECAMERAVIEWMODERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean UndoView()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_UNDOVIEWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RedoView()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_REDOVIEWRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleLockViewMode()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLELOCKVIEWMODERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleLockTools()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLELOCKTOOLSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleLockTime()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLELOCKTIMERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetWindowMode(int mode)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETWINDOWMODERPC);
        rpc.SetWindowMode(mode);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ToggleBoundingBoxMode()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_TOGGLEBOUNDINGBOXMODERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean EnableTool(int tool, boolean enabled)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ENABLETOOLRPC);
        rpc.SetToolId(tool);
        rpc.SetBoolFlag(enabled);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CopyViewToWindow(int from, int to)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_COPYVIEWTOWINDOWRPC);
        rpc.SetWindowLayout(from);
        rpc.SetWindowId(to);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CopyLightingToWindow(int from, int to)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_COPYLIGHTINGTOWINDOWRPC);
        rpc.SetWindowLayout(from);
        rpc.SetWindowId(to);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CopyAnnotationsToWindow(int from, int to)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_COPYANNOTATIONSTOWINDOWRPC);
        rpc.SetWindowLayout(from);
        rpc.SetWindowId(to);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean CopyPlotsToWindow(int from, int to)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_COPYPLOTSTOWINDOWRPC);
        rpc.SetWindowLayout(from);
        rpc.SetWindowId(to);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetAnnotationAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETANNOTATIONATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultAnnotationAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTANNOTATIONATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetAnnotationAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETANNOTATIONATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean AddAnnotationObject(int annotType)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_ADDANNOTATIONOBJECTRPC);
        rpc.SetIntArg1(annotType);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean HideActiveAnnotationObjects()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_HIDEACTIVEANNOTATIONOBJECTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DeleteActiveAnnotationObjects()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DELETEACTIVEANNOTATIONOBJECTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean RaiseActiveAnnotationObjects()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RAISEACTIVEANNOTATIONOBJECTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean LowerActiveAnnotationObjects()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_LOWERACTIVEANNOTATIONOBJECTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetAnnotationObjectOptions()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETANNOTATIONOBJECTOPTIONSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultAnnotationObjectList()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTANNOTATIONOBJECTLISTRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetAnnotationObjectList()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETANNOTATIONOBJECTLISTRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetLightList()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETLIGHTLISTRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultLightList()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTLIGHTLISTRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetLightList()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETLIGHTLISTRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetAnimationAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETANIMATIONATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetAppearanceAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETAPPEARANCERPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ProcessExpressions()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_PROCESSEXPRESSIONSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearPickPoints()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARPICKPOINTSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ClearReferenceLines()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_CLEARREFLINESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetViewExtents(int t)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETVIEWEXTENTSTYPERPC);
        rpc.SetWindowLayout(t);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DatabaseQuery(String queryName, Vector vars)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DATABASEQUERYRPC);
        rpc.SetQueryName(queryName);
        rpc.SetQueryVariables(vars);
        rpc.SetIntArg1(0);
        rpc.SetIntArg2(0);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean DatabaseQuery(String queryName, Vector vars,
                   int arg1, int arg2)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_DATABASEQUERYRPC);
        rpc.SetQueryName(queryName);
        rpc.SetQueryVariables(vars);
        rpc.SetIntArg1(arg1);
        rpc.SetIntArg2(arg2);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean PointQuery(String queryName, double[] pt, Vector vars)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_POINTQUERYRPC);
        rpc.SetQueryName(queryName);
        rpc.SetQueryPoint1(pt);
        rpc.SetQueryVariables(vars);
        rpc.SetIntArg1(-1);
        rpc.SetIntArg2(-1);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean PointQuery(String queryName, double[] pt, Vector vars,
                              int arg1, int arg2)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_POINTQUERYRPC);
        rpc.SetQueryName(queryName);
        rpc.SetQueryPoint1(pt);
        rpc.SetQueryVariables(vars);
        rpc.SetIntArg1(arg1);
        rpc.SetIntArg2(arg2);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean LineQuery(String queryName, double[] pt1, double[] pt2,
        Vector vars)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_LINEQUERYRPC);
        rpc.SetQueryName(queryName);
        rpc.SetQueryPoint1(pt1);
        rpc.SetQueryPoint2(pt2);
        rpc.SetQueryVariables(vars);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean Pick(int x, int y, Vector vars)
    {
        double[] pt = new double[3];
        pt[0] = (double)x;
        pt[1] = (double)y;
        pt[2] = 0.;
        return PointQuery("ScreenZonePick", pt, vars);
    }

    public boolean Lineout(double x0, double y0, double x1, double y1, Vector vars)
    {
        double[] pt1 = new double[3];
        pt1[0] = x0;
        pt1[1] = y0;
        pt1[2] = 0.;
        double[] pt2 = new double[3];
        pt2[0] = x1;
        pt2[1] = y1;
        pt2[2] = 0.;
        return LineQuery("Lineout", pt1, pt2, vars);
    }

    public boolean SetMaterialAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETMATERIALATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultMaterialAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTMATERIALATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetMaterialAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETMATERIALATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetKeyframeAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETKEYFRAMEATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetGlobalLineoutAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETGLOBALLINEOUTATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetPickAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETPICKATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultPickAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTPICKATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetPickAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETPICKATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetPickLetter()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETPICKLETTERRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetQueryOverTimeAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETQUERYOVERTIMEATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultQueryOverTimeAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTQUERYOVERTIMEATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetQueryOverTimeAttributes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETQUERYOVERTIMEATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetInteractorAttriubes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETINTERACTORATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetDefaultInteractorAttriubes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETDEFAULTINTERACTORATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResetQueryInteractorAttriubes()
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESETINTERACTORATTRIBUTESRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ShowToolbars(boolean forAllWindows)
    {
        if(forAllWindows)
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SHOWTOOLBARSFORALLWINDOWSRPC);
        else
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SHOWTOOLBARSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean HideToolbars(boolean forAllWindows)
    {
        if(forAllWindows)
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_HIDETOOLBARSFORALLWINDOWSRPC);
        else
            rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_HIDETOOLBARSRPC);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean OpenClient(String clientName, String program, Vector args)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_OPENCLIENTRPC);
        rpc.SetDatabase(clientName);
        rpc.SetProgramHost(program);
        rpc.SetProgramOptions(args);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean SetTryHarderCyclesTimes(int flag)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_SETTRYHARDERCYCLESTIMESRPC);
        rpc.SetIntArg1(flag);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean ResizeWindow(int win, int w, int h)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_RESIZEWINDOWRPC);
        rpc.SetWindowId(win);
        rpc.SetIntArg1(w);
        rpc.SetIntArg2(h);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean MoveWindow(int win, int x, int y)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_MOVEWINDOWRPC);
        rpc.SetWindowId(win);
        rpc.SetIntArg1(x);
        rpc.SetIntArg2(y);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public boolean MoveAndResizeWindow(int win, int x, int y, int w, int h)
    {
        rpc.SetRPCType(ViewerRPC.VIEWERRPCTYPE_MOVEANDRESIZEWINDOWRPC);
        rpc.SetWindowId(win);
        rpc.SetIntArg1(x);
        rpc.SetIntArg2(y);
        rpc.SetIntArg3(w);
        rpc.SetWindowLayout(h);
        rpc.Notify();
        return synchronous ? Synchronize() : true;
    }

    public synchronized boolean Synchronize()
    {
        // Clear any error in the message observer.
        messageObserver.ClearError();

        // Tell the syncObserver to notify this thread when the new syncCount
        // gets returned from the viewer.
        ++syncCount;
        syncNotifier.NotifyThreadOnValue(Thread.currentThread(), syncCount);
        syncNotifier.SetUpdate(false);

        // Send a new synchronization tag
        syncAtts.SetSyncTag(syncCount);
        syncAtts.Notify();
        syncAtts.SetSyncTag(-1);

        // Wait until the syncNotifier notifies us that it's okay to proceed.
        boolean noErrors = true;
        try
        {
            synchronized(Thread.currentThread())
            { 
                PrintMessage("Waiting for viewer sync.");
                Thread.currentThread().wait();
            }
        }
        catch(java.lang.InterruptedException e)
        {
            noErrors = false;
        }

        // Make it return true if there are no errors.
        return noErrors && !messageObserver.GetErrorFlag();
    }

    public String GetLastError()
    {
        return messageObserver.GetLastError();
    }

//
// Methods to get the state objects
//
    public MessageAttributes GetMessageAttributes() { return messageAtts; }
    public AnnotationAttributes GetAnnotationAttributes() { return annotationAtts; }
    public GlobalAttributes GetGlobalAttributes() { return globalAtts;}
    public DatabaseCorrelationList GetDatabaseCorrelationList() { return correlationList; }
    public PlotList GetPlotList() { return plotList; }
    public PluginManagerAttributes GetPluginAtts() { return pluginAtts;}
    public HostProfileList GetHostProfiles() { return hostProfiles; }
    public AppearanceAttributes GetAppearanceAttributes() { return appearanceAtts; }
    public ViewCurveAttributes GetViewCurve() { return viewCurve; }
    public View2DAttributes GetView2D() { return view2D; }
    public View3DAttributes GetView3D() { return view3D; }
    public ColorTableAttributes GetColorTableAttributes() { return colorTableAtts; }
    public SaveWindowAttributes GetSaveWindowAttributes() { return saveAtts; }
    public StatusAttributes GetStatusAttributes() { return statusAtts; }
    public EngineList GetEngineList() { return engineList; }
    public LightList GetLightList() { return lightList; }
    public AnimationAttributes GetAnimationAttributes() { return animationAtts; }
    public ExpressionList GetExpressionList() { return expressionList; }
    public PickAttributes GetPickAttributes() { return pickAtts; }
    public WindowInformation GetWindowInformation() { return windowInfo; }
    public RenderingAttributes GetRenderingAttributes() { return renderAtts; }
    public QueryList GetQueryList() { return queryList; }
    public QueryAttributes GetQueryAttributes() { return queryAtts; }
    public MaterialAttributes GetMaterialAttributes() { return materialAtts; }
    public GlobalLineoutAttributes GetGlobalLineoutAttributes() { return globalLineoutAtts; }
    public AnnotationObjectList GetAnnotationObjectList() { return annotationObjectList; }
    public QueryOverTimeAttributes GetQueryOverTimeAttributes() { return queryOverTimeAtts; }
    public InteractorAttributes GetInteractorAttributes() { return interactorAtts; }
    public ProcessAttributes GetProcessAttributes() { return processAtts; }
    public DBPluginInfoAttributes GetPluginInfoAttributes() { return dbPluginInfoAtts; }
    public ExportDBAttributes GetExportDBAttributes() { return exportDBAtts; }
    public ClientMethod GetClientMethod() { return clientMethod; }
    public ClientInformation GetClientInformation() { return clientInformation; }
    public final ClientInformationList GetClientInformationList() { return clientInformationList; }
    public MovieAttributes GetMovieAttributes() { return movieAtts; }
    public MeshManagementAttributes GetMeshManagementAttributes() { return meshManagementAtts; }
    public ConstructDDFAttributes GetDDFAttributes() { return constructDDFAtts; }
    public PlotInfoAttributes GetPlotInfoAttributes() { return plotInfoAtts; }

    public int GetPlotIndex(String plotName)
    {
        return plotPlugins.IndexFromName(plotName);
    }
 
    public String GetPlotName(int index) throws ArrayIndexOutOfBoundsException
    { 
        return plotPlugins.GetPluginName(index);
    }

    public String GetPlotVersion(int index) throws ArrayIndexOutOfBoundsException
    { 
        return plotPlugins.GetPluginVersion(index);
    }

    public Plugin GetPlotAttributes(int index) throws ArrayIndexOutOfBoundsException
    {
        return plotPlugins.GetPluginAttributes(index);
    }

    public Plugin GetPlotAttributes(String plotName) throws ArrayIndexOutOfBoundsException
    {
        int index = plotPlugins.IndexFromName(plotName);
        return plotPlugins.GetPluginAttributes(index);
    }

    public int GetNumPlotPlugins()
    {
        return plotPlugins.GetNumPlugins();
    }

    public int GetOperatorIndex(String plotName)
    {
        return operatorPlugins.IndexFromName(plotName);
    }

    public String GetOperatorName(int index) throws ArrayIndexOutOfBoundsException
    { 
        return operatorPlugins.GetPluginName(index);
    }

    public String GetOperatorVersion(int index) throws ArrayIndexOutOfBoundsException
    { 
        return operatorPlugins.GetPluginVersion(index);
    }
    
    public Plugin GetOperatorAttributes(int index) throws ArrayIndexOutOfBoundsException
    {
        return operatorPlugins.GetPluginAttributes(index);
    }

    public Plugin GetOperatorAttributes(String operatorName) throws ArrayIndexOutOfBoundsException
    {
        int index = operatorPlugins.IndexFromName(operatorName);
        return operatorPlugins.GetPluginAttributes(index);
    }

    public int GetNumOperatorPlugins()
    {
        return operatorPlugins.GetNumPlugins();
    }

//
// SimpleObserver interface
//
    public void Update(AttributeSubject s)
    {
        if(s == pluginAtts)
        {
            LoadPlugins();
        }
    }

    public void SetUpdate(boolean val) { doUpdate = val; }
    public boolean GetUpdate() { return doUpdate; }

    public boolean PluginsLoaded()
    {
        return pluginsLoaded;
    }

    private void LoadPlugins()
    {
        if(!pluginsLoaded)
        {
            // Try loading the plot plugins. If they can all be loaded,
            // add them to xfer.
            plotPlugins.LoadPlugins(pluginAtts);
            for(int i = 0; i < plotPlugins.GetNumPlugins(); ++i)
            {
                 Plugin p = plotPlugins.GetPluginAttributes(i);
                 if(p != null)
                     xfer.Add((AttributeSubject)p);
                 else
                     xfer.AddDummy();
            }

            // Try loading the operator plugins. If they can all be loaded,
            // add them to xfer.
            operatorPlugins.LoadPlugins(pluginAtts);
            for(int i = 0; i < operatorPlugins.GetNumPlugins(); ++i)
            {
                 Plugin p = operatorPlugins.GetPluginAttributes(i);
                 if(p != null)
                     xfer.Add((AttributeSubject)p);
                 else
                     xfer.AddDummy();
            }

            pluginsLoaded = true;
        }
    }
//
// Data members
//
    private RemoteProcess        viewer;
    private Xfer                 xfer;
    private PluginManager        plotPlugins;
    private PluginManager        operatorPlugins;
    private int                  syncCount;
    private boolean              synchronous;
    private MessageObserver      messageObserver;
    private SyncNotifier         syncNotifier;
    private boolean              pluginsLoaded;
    private boolean              doUpdate;
    private boolean              verbose;
    private String               dataPath;

//
// State objects
//
    private ViewerRPC                rpc;
    private PostponedAction          postponedAction;
    private SyncAttributes           syncAtts;
    private AppearanceAttributes     appearanceAtts;
    private PluginManagerAttributes  pluginAtts;
    private GlobalAttributes         globalAtts;
    private DatabaseCorrelationList  correlationList;
    private PlotList                 plotList;
    private HostProfileList          hostProfiles;
    private MessageAttributes        messageAtts;
    private SaveWindowAttributes     saveAtts;
    private StatusAttributes         statusAtts;
    private EngineList               engineList;
    private ColorTableAttributes     colorTableAtts;
    private ExpressionList           expressionList;
    private AnnotationAttributes     annotationAtts;
    private SILRestrictionAttributes silRestrictionAtts;
    private ViewCurveAttributes      viewCurve;
    private View2DAttributes         view2D;
    private View3DAttributes         view3D;
    private LightList                lightList;
    private AnimationAttributes      animationAtts;
    private PickAttributes           pickAtts;
    private PrinterAttributes        printerAtts;
    private WindowInformation        windowInfo;
    private RenderingAttributes      renderAtts;
    private KeyframeAttributes       keyframeAtts;
    private QueryList                queryList;
    private QueryAttributes          queryAtts;
    private MaterialAttributes       materialAtts;
    private GlobalLineoutAttributes  globalLineoutAtts;
    private AnnotationObjectList     annotationObjectList;
    private QueryOverTimeAttributes  queryOverTimeAtts;
    private InteractorAttributes     interactorAtts;
    private SILRestrictionAttributes silAtts;
    private ProcessAttributes        processAtts;
    private DBPluginInfoAttributes   dbPluginInfoAtts;
    private ExportDBAttributes       exportDBAtts;
    private ConstructDDFAttributes   constructDDFAtts;
    private ClientMethod             clientMethod;
    private ClientInformation        clientInformation;
    private ClientInformationList    clientInformationList;
    private MovieAttributes          movieAtts;
    private MeshManagementAttributes meshManagementAtts;
    private ViewerRPC                logRPC;
    private PlotInfoAttributes       plotInfoAtts;
}
