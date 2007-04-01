// ****************************************************************************
//                               TimeVaryingExodusPluginInfo.h
// ****************************************************************************

#ifndef TIMEVARYINGEXODUS_PLUGIN_INFO_H
#define TIMEVARYINGEXODUS_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: TimeVaryingExodusDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the TimeVaryingExodus plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Tue Feb 22 14:48:28 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class TimeVaryingExodusGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class TimeVaryingExodusCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual TimeVaryingExodusGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class TimeVaryingExodusMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual TimeVaryingExodusCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class TimeVaryingExodusEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual TimeVaryingExodusCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
