// ****************************************************************************
//                               ANSYSPluginInfo.h
// ****************************************************************************

#ifndef ANSYS_PLUGIN_INFO_H
#define ANSYS_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: ANSYSDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the ANSYS plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Wed Jul 6 22:10:40 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class ANSYSGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class ANSYSCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual ANSYSGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class ANSYSMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual ANSYSCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class ANSYSEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual ANSYSCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
