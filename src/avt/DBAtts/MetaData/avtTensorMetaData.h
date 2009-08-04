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

#ifndef AVTTENSORMETADATA_H
#define AVTTENSORMETADATA_H
#include <dbatts_exports.h>
#include <string>
#include <avtTypes.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: avtTensorMetaData
//
// Purpose:
//    Contains tensor metadata attributes
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class DBATTS_API avtTensorMetaData : public AttributeSubject
{
public:
    avtTensorMetaData();
    avtTensorMetaData(const avtTensorMetaData &obj);
    virtual ~avtTensorMetaData();

    virtual avtTensorMetaData& operator = (const avtTensorMetaData &obj);
    virtual bool operator == (const avtTensorMetaData &obj) const;
    virtual bool operator != (const avtTensorMetaData &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // User-defined methods
    avtTensorMetaData(std::string, std::string, avtCentering, int);
    void Print(ostream &, int = 0) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_name = 0,
        ID_originalName,
        ID_validVariable,
        ID_meshName,
        ID_centering,
        ID_hasUnits,
        ID_units,
        ID_dim,
        ID_hideFromGUI
    };

public:
    std::string  name;
    std::string  originalName;
    bool         validVariable;
    std::string  meshName;
    avtCentering centering;
    bool         hasUnits;
    std::string  units;
    int          dim;
    bool         hideFromGUI;

private:
    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
