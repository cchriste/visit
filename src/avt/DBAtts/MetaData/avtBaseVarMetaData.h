/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400124
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
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
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef AVTBASEVARMETADATA_H
#define AVTBASEVARMETADATA_H
#include <dbatts_exports.h>
#include <string>
#include <AttributeSubject.h>

#include <vector>

// ****************************************************************************
// Class: avtBaseVarMetaData
//
// Purpose:
//    Contains metadata attributes associated with all mesh variables
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class DBATTS_API avtBaseVarMetaData : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    avtBaseVarMetaData();
    avtBaseVarMetaData(const avtBaseVarMetaData &obj);
protected:
    // These constructors are for objects derived from this class
    avtBaseVarMetaData(private_tmfs_t tmfs);
    avtBaseVarMetaData(const avtBaseVarMetaData &obj, private_tmfs_t tmfs);
public:
    virtual ~avtBaseVarMetaData();

    virtual avtBaseVarMetaData& operator = (const avtBaseVarMetaData &obj);
    virtual bool operator == (const avtBaseVarMetaData &obj) const;
    virtual bool operator != (const avtBaseVarMetaData &obj) const;
private:
    void Init();
    void Copy(const avtBaseVarMetaData &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // User-defined methods
    avtBaseVarMetaData(private_tmfs_t, std::string, std::string);
    void Print(ostream &, int = 0) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_name = 0,
        ID_originalName,
        ID_meshName,
        ID_validVariable,
        ID_hideFromGUI,
        ID__LAST
    };

public:
    std::string name;
    std::string originalName;
    std::string meshName;
    bool        validVariable;
    bool        hideFromGUI;

private:
    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define AVTBASEVARMETADATA_TMFS "sssbb"

#endif
