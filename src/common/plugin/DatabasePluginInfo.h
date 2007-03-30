// ************************************************************************* //
//                              DatabasePluginInfo.h                             //
// ************************************************************************* //

#ifndef DATABASE_PLUGIN_INFO_H
#define DATABASE_PLUGIN_INFO_H
#include <plugin_exports.h>
#include <stdio.h>

#include <string>
#include <vector>

enum DatabaseType
{
    DB_TYPE_STSD,
    DB_TYPE_STMD,
    DB_TYPE_MTSD,
    DB_TYPE_MTMD,
    DB_TYPE_CUSTOM
};

// Forward declarations.
class avtDatabase;

// ****************************************************************************
//  Class: *DatabasePluginInfo
//
//  Purpose:
//    Classes that provide all the information about the database plugin.
//    Portions are separated into pieces relevant to the appropriate
//    components of VisIt.
//
//  Programmer: Jeremy Meredith
//  Creation:   August 21, 2002
//
//  Modifications:
//
// ****************************************************************************

class PLUGIN_API GeneralDatabasePluginInfo
{
  public:
    virtual char *GetName() const = 0;
    virtual char *GetVersion() const = 0;
    virtual char *GetID() const = 0;
};

class PLUGIN_API CommonDatabasePluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
    virtual DatabaseType              GetDatabaseType() = 0;
    virtual std::vector<std::string>  GetDefaultExtensions() = 0;
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock) = 0;
};

class PLUGIN_API MDServerDatabasePluginInfo : public virtual CommonDatabasePluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy() = 0;
};

class PLUGIN_API EngineDatabasePluginInfo : public virtual CommonDatabasePluginInfo
{
  public:
    // this makes compilers happy... remove if we ever have functions here
    virtual void dummy() = 0;
};

#endif
