// ****************************************************************************
//                               ShapefilePluginInfo.h
// ****************************************************************************

#ifndef SHAPEFILE_PLUGIN_INFO_H
#define SHAPEFILE_PLUGIN_INFO_H
#include <DatabasePluginInfo.h>
#include <database_plugin_exports.h>

class avtDatabase;
class avtDatabaseWriter;

// ****************************************************************************
//  Class: ShapefileDatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the Shapefile plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: whitlocb -- generated by xml2info
//  Creation:   Mon Mar 28 03:13:03 PDT 2005
//
//  Modifications:
//
// ****************************************************************************

class ShapefileGeneralPluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const;
    virtual char *GetVersion() const;
    virtual char *GetID() const;
    virtual bool  EnabledByDefault() const;
    virtual bool  HasWriter() const;
};

class ShapefileCommonPluginInfo : public virtual CommonDatabasePluginInfo, public virtual ShapefileGeneralPluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType();
    virtual std::vector<std::string>  GetDefaultExtensions();
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock);
};

class ShapefileMDServerPluginInfo : public virtual MDServerDatabasePluginInfo, public virtual ShapefileCommonPluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy();
};

class ShapefileEnginePluginInfo : public virtual EngineDatabasePluginInfo, public virtual ShapefileCommonPluginInfo
{
  public:
    virtual avtDatabaseWriter        *GetWriter(void);
};

#endif
