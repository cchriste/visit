// ****************************************************************************
//                               Boxlib3DPluginInfo.h
// ****************************************************************************

#ifndef BOXLIB3D_PLUGIN_INFO_H
#define BOXLIB3D_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: Boxlib3DDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the Boxlib3D plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Thu Nov 6 10:33:13 PDT 2003
//
//  Modifications:
//
// ****************************************************************************

class Boxlib3DGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
};

class Boxlib3DCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual Boxlib3DGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
    virtual avtDatabaseWriter        *GetWriter(void);
};

class Boxlib3DMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual Boxlib3DCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class Boxlib3DEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual Boxlib3DCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

#endif
