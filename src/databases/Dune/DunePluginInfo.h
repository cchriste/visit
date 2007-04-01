// ****************************************************************************
//                               DunePluginInfo.h
// ****************************************************************************

#ifndef DUNE_PLUGIN_INFO_H
#define DUNE_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: DuneDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the Dune plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Tue Feb 22 14:36:58 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class DuneGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class DuneCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual DuneGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class DuneMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual DuneCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class DuneEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual DuneCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
