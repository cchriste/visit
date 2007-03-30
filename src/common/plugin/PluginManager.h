// ************************************************************************* //
//                             PluginManager.h                               //
// ************************************************************************* //

#ifndef PLUGIN_MANAGER_H
#define PLUGIN_MANAGER_H
#include <plugin_exports.h>
#include <string>
#include <vector>
#include <map>
#include <utility>

// ****************************************************************************
//  Class: PluginManager
//
//  Purpose:
//    The plugin manager.  It provides an abstraction for all plugin
//    managers.  The information is broken up into several classes since
//    portions of it are only relevant to particular components within visit.
//    There is the general information which all the components are interested
//    in, then portions for the gui, viewer, cli, engine, and mdserver.
//
//  Programmer: Jeremy Meredith
//  Creation:   August 20, 2002
//
//  Modifications:
//    Jeremy Meredith, Fri Feb 28 12:28:50 PST 2003
//    Renamed some methods and data members to make their function and
//    usage more correct and obvious.  Added support for loading plugins
//    on demand.  Made PluginLoaded be private and added PluginAvailable,
//    which can attempt to load a plugin on demand before checking to see
//    if it is loaded.
//
// ****************************************************************************

class PLUGIN_API PluginManager
{
  public:
    enum PluginCategory
    {
        no_category,
        GUI,
        Viewer,
        Engine,
        MDServer,
        Scripting
    };
  public:
    virtual ~PluginManager();

    static void                     Initialize(const std::string &managerName,
                                               const PluginCategory,
                                               bool=false);

    void                            DisablePlugin(const std::string&);
    void                            EnablePlugin(const std::string&);

    virtual void                    LoadPluginsNow();
    virtual void                    LoadPluginsOnDemand();
    virtual void                    ReloadPlugins();
    virtual void                    UnloadPlugins();

    bool                            PluginExists(const std::string&);
    bool                            PluginAvailable(const std::string&);

    std::string                     GetPluginName(const std::string&);
    std::string                     GetPluginVersion(const std::string&);

    int                             GetNAllPlugins() const;
    std::string                     GetAllID(const int) const;
    int                             GetAllIndex(const std::string &) const;
    int                             GetAllIndexFromName(const std::string &) const;

    int                             GetNEnabledPlugins() const;
    std::string                     GetEnabledID(const int) const;

  protected:
                                    PluginManager(const std::string&);
    void                            ReadPluginInfo();
    void                            SetPluginDir();
    void                            ReadPluginDir(std::vector<
                                                  std::pair<std::string,
                                                            std::string> > &);
    void                            GetPluginList(std::vector<
                                                  std::pair<std::string,
                                                             std::string> >&);
    bool                            IsGeneralPlugin(const std::string &) const;

    bool                            PluginLoaded(const std::string&);
    void                            PluginOpen(const std::string &pluginFile);
    void                           *PluginSymbol(const std::string &symbol);
    char                           *PluginError() const;
    void                            PluginClose();

    virtual void                    LoadSinglePlugin(int i);

    virtual bool                    LoadGeneralPluginInfo()    = 0;
    virtual void                    LoadGUIPluginInfo()        { }
    virtual void                    LoadViewerPluginInfo()     { }
    virtual void                    LoadMDServerPluginInfo()   { }
    virtual void                    LoadEnginePluginInfo()     { }
    virtual void                    LoadScriptingPluginInfo()  { }

    virtual void                    FreeCommonPluginInfo()     = 0;
    virtual void                    FreeGUIPluginInfo()        { }
    virtual void                    FreeViewerPluginInfo()     { }
    virtual void                    FreeMDServerPluginInfo()   { }
    virtual void                    FreeEnginePluginInfo()     { }
    virtual void                    FreeScriptingPluginInfo()  { }

    std::vector<std::string>                pluginDirs;
    std::string                             openPlugin;
    void                                   *handle;
    char                                   *pluginError;
    int                                     category;
    bool                                    parallel;
    std::string                             managerName;
    bool                                    loadOnDemand;

    // arrays containing all plugins
    std::vector<std::string>                ids;
    std::vector<std::string>                names;
    std::vector<std::string>                versions;
    std::vector<std::string>                libfiles;
    std::vector<bool>                       enabled;

    // maps from id->allindex and id->loadedindex
    std::map<std::string, int>              allindexmap;
    std::map<std::string, int>              loadedindexmap;

    // arrays containing enabled plugins
    std::vector<void*>                      loadedhandles;
    std::vector<std::string>                loadedids;
};

#endif
