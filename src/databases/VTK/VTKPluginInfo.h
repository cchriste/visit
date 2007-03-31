// ****************************************************************************
//                               VTKPluginInfo.h
// ****************************************************************************

#ifndef VTK_PLUGIN_INFO_H
#define VTK_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: VTKDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the VTK plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Fri Sep 12 15:03:44 PST 2003
//
//  Modifications:
//
// ****************************************************************************

class VTKGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class VTKCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual VTKGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
    virtual avtDatabaseWriter        *GetWriter(void);
};

class VTKMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual VTKCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class VTKEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual VTKCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

#endif
