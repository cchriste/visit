// ****************************************************************************
//                               EnSightPluginInfo.h
// ****************************************************************************

#ifndef ENSIGHT_PLUGIN_INFO_H
#define ENSIGHT_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: EnSightDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the EnSight plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Tue Feb 22 14:37:19 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class EnSightGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class EnSightCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual EnSightGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class EnSightMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual EnSightCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class EnSightEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual EnSightCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
