// ****************************************************************************
//                               KullLitePluginInfo.h
// ****************************************************************************

#ifndef KULLLITE_PLUGIN_INFO_H
#define KULLLITE_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: KullLiteDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the KullLite plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Mon Mar 22 09:18:32 PDT 2004
//
//  Modifications:
//
// ****************************************************************************

class KullLiteGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
};

class KullLiteCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual KullLiteGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
    virtual avtDatabaseWriter        *GetWriter(void);
};

class KullLiteMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual KullLiteCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class KullLiteEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual KullLiteCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

#endif
