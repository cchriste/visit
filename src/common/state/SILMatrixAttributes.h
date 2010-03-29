/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

#ifndef SILMATRIXATTRIBUTES_H
#define SILMATRIXATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: SILMatrixAttributes
//
// Purpose:
//    This class contain the information needed to represent a SIL Matrix.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SILMatrixAttributes : public AttributeSubject
{
public:
    SILMatrixAttributes();
    SILMatrixAttributes(const SILMatrixAttributes &obj);
    virtual ~SILMatrixAttributes();

    virtual SILMatrixAttributes& operator = (const SILMatrixAttributes &obj);
    virtual bool operator == (const SILMatrixAttributes &obj) const;
    virtual bool operator != (const SILMatrixAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSet1();
    void SelectCategory1();
    void SelectSet2();
    void SelectCategory2();

    // Property setting methods
    void SetSet1(const intVector &set1_);
    void SetCategory1(const std::string &category1_);
    void SetRole1(int role1_);
    void SetSet2(const intVector &set2_);
    void SetCategory2(const std::string &category2_);
    void SetRole2(int role2_);

    // Property getting methods
    const intVector   &GetSet1() const;
          intVector   &GetSet1();
    const std::string &GetCategory1() const;
          std::string &GetCategory1();
    int               GetRole1() const;
    const intVector   &GetSet2() const;
          intVector   &GetSet2();
    const std::string &GetCategory2() const;
          std::string &GetCategory2();
    int               GetRole2() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_set1 = 0,
        ID_category1,
        ID_role1,
        ID_set2,
        ID_category2,
        ID_role2
    };

private:
    intVector   set1;
    std::string category1;
    int         role1;
    intVector   set2;
    std::string category2;
    int         role2;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
