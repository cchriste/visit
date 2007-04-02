/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            DatabasePluginInfo.h                           //
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
class avtDatabaseWriter;
class DBOptionsAttributes;

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
//    Hank Childs, Wed Sep 10 07:05:54 PDT 2003
//    Added DatabaseWriter.
//
//    Jeremy Meredith, Wed Nov  5 10:28:29 PST 2003
//    Added ability to disable plugins by default.
//
//    Hank Childs, Thu Feb 19 10:01:47 PST 2004
//    Added GetFilenames.  Made GetDefaultExtensions not be pure virtual.
//
//    Jeremy Meredith, Tue Feb 22 18:36:54 PST 2005
//    Moved GetWriter to the engine so the mdserver doesn't need it.
//    Added the general plugin info method HasWriter so the mdserver
//    can still check if it is supported by the given plugin.
//
//    Hank Childs, Tue Mar 22 16:06:15 PST 2005
//    Make destructor virtual.
//
//    Hank Childs, Mon May 23 16:31:36 PDT 2005
//    Add DBOptions.
//
// ****************************************************************************

class PLUGIN_API GeneralDatabasePluginInfo
{
  public:
    virtual ~GeneralDatabasePluginInfo() {;};
    virtual char *GetName() const = 0;
    virtual char *GetVersion() const = 0;
    virtual char *GetID() const = 0;
    virtual bool  EnabledByDefault() const { return true; }
    virtual bool  HasWriter() const { return false; }
};

class PLUGIN_API CommonDatabasePluginInfo : public virtual GeneralDatabasePluginInfo
{
  public:
                                      CommonDatabasePluginInfo();
    virtual                          ~CommonDatabasePluginInfo();

    virtual DatabaseType              GetDatabaseType() = 0;
    virtual std::vector<std::string>  GetDefaultExtensions()
                                   { std::vector<std::string> rv; return rv; };
    virtual std::vector<std::string>  GetFilenames()
                                   { std::vector<std::string> rv; return rv; };
    virtual avtDatabase              *SetupDatabase(const char * const *list,
                                                    int nList, int nBlock) = 0;

    virtual DBOptionsAttributes      *GetReadOptions(void) const;
    virtual DBOptionsAttributes      *GetWriteOptions(void) const;
    virtual void                      SetReadOptions(DBOptionsAttributes *);
    virtual void                      SetWriteOptions(DBOptionsAttributes *);

  protected:
    DBOptionsAttributes              *readOptions;
    DBOptionsAttributes              *writeOptions;
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
    virtual avtDatabaseWriter        *GetWriter(void) { return NULL; };
};

#endif
