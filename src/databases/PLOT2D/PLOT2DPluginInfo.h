// ****************************************************************************
//                               PLOT2DPluginInfo.h
// ****************************************************************************

#ifndef PLOT2D_PLUGIN_INFO_H
#define PLOT2D_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: PLOT2DDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the PLOT2D plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: meredith -- generated by xml2info
//  Creation:   Tue Feb 22 14:41:37 PST 2005
//
//  Modifications:
//
// ****************************************************************************

class PLOT2DGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class PLOT2DCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual PLOT2DGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class PLOT2DMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual PLOT2DCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class PLOT2DEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual PLOT2DCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
