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
//                           PlotPluginManager.h                             //
// ************************************************************************* //

#ifndef PLOT_PLUGIN_MANAGER_H
#define PLOT_PLUGIN_MANAGER_H
#include <plugin_exports.h>
#include <PluginManager.h>
#include <string>
#include <vector>

class GeneralPlotPluginInfo;
class CommonPlotPluginInfo;
class GUIPlotPluginInfo;
class ViewerPlotPluginInfo;
class EnginePlotPluginInfo;
class ScriptingPlotPluginInfo;

// ****************************************************************************
//  Class: PlotPluginManager
//
//  Purpose:
//    The plot plugin manager.  It reads in the operator
//    plugins and provides information about the plot plugins.
//
//  Programmer: Eric Brugger
//  Creation:   March 7, 2001
//
//  Modifications:
//    Jeremy Meredith, Fri May 11 14:09:33 PDT 2001
//    Added a PluginCategory which must be passed to Initialize.
//    Added an array of handles for the real shared libraries.
//
//    Jeremy Meredith, Thu Jul 26 09:32:55 PDT 2001
//    Renamed plugin info classes.
//
//    Jeremy Meredith, Thu Sep  6 16:23:17 PDT 2001
//    Reworked this to support more advanced plugin management, and 
//    support for runtime plugin loading/unloading.
//
//    Brad Whitlock, Thur Jul 12 17:44:54 PST 2001
//    Added support for scripting plugins.
//
//    Jeremy Meredith, Tue Jan 15 16:12:47 PST 2002
//    Added runtime parallel config.
//
//    Brad Whitlock, Fri Mar 29 09:21:27 PDT 2002
//    Ported to windows.
//
//    Jeremy Meredith, Fri Jul  5 17:36:23 PDT 2002
//    Made it work on plugins from multiple directories.
//
//    Jeremy Meredith, Tue Jul  9 17:42:07 PDT 2002
//    Made LoadCommonPluginInfo return a success flag.
//
//    Jeremy Meredith, Tue Aug 20 16:57:40 PDT 2002
//    Refactored most stuff into common PluginManager base class.
//
//    Jeremy Meredith, Fri Feb 28 12:26:45 PST 2003
//    Renamed LoadCommonPluginInfo to LoadGeneralPluginInfo to clarify its
//    purpose (which really is unrelated to FreeCommonPluginInfo and 
//    GetCommonPluginInfo).
//
//    Brad Whitlock, Wed Nov 22 16:30:48 PST 2006
//    Added const char * argument to Initialize function.
//
//    Brad Whitlock, Tue Jun 24 10:48:11 PDT 2008
//    Removed singleton characteristics.
//
//    Brad Whitlock, Wed Jun 17 10:24:46 PDT 2009
//    I added arguments to Initialize.
//
// ****************************************************************************

class PLUGIN_API PlotPluginManager : public PluginManager
{
  public:
                                    PlotPluginManager();
    virtual                        ~PlotPluginManager();

    void                            Initialize(const PluginCategory,
                                               bool parallel=false,
                                               const char *plugindir = 0,
                                               bool readInfo = true,
                                               PluginBroadcaster *broadcaster=0);

    virtual void                    ReloadPlugins();

    CommonPlotPluginInfo           *GetCommonPluginInfo(const std::string&);
    GUIPlotPluginInfo              *GetGUIPluginInfo(const std::string&);
    ViewerPlotPluginInfo           *GetViewerPluginInfo(const std::string&);
    EnginePlotPluginInfo           *GetEnginePluginInfo(const std::string&);
    ScriptingPlotPluginInfo        *GetScriptingPluginInfo(const std::string&);

  private:
    virtual bool                    LoadGeneralPluginInfo();
    virtual void                    LoadGUIPluginInfo();
    virtual void                    LoadViewerPluginInfo();
    virtual void                    LoadEnginePluginInfo();
    virtual void                    LoadScriptingPluginInfo();

    virtual void                    FreeCommonPluginInfo();
    virtual void                    FreeGUIPluginInfo();
    virtual void                    FreeViewerPluginInfo();
    virtual void                    FreeEnginePluginInfo();
    virtual void                    FreeScriptingPluginInfo();

    // arrays containing enabled plugins
    std::vector<CommonPlotPluginInfo*>      commonPluginInfo;
    std::vector<GUIPlotPluginInfo*>         guiPluginInfo;
    std::vector<ViewerPlotPluginInfo*>      viewerPluginInfo;
    std::vector<EnginePlotPluginInfo*>      enginePluginInfo;
    std::vector<ScriptingPlotPluginInfo*>   scriptingPluginInfo;
};

#endif
