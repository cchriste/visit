// ****************************************************************************
//                               MFIXPluginInfo.h
// ****************************************************************************

#ifndef MFIX_PLUGIN_INFO_H
#define MFIX_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: MFIXDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the MFIX plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: bdotson -- generated by xml2info
//  Creation:   Fri May 26 08:59:22 PDT 2006
//
//  Modifications:
//
// ****************************************************************************

class MFIXGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class MFIXCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual MFIXGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class MFIXMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual MFIXCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class MFIXEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual MFIXCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
