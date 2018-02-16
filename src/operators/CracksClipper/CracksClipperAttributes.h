/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
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

#ifndef CRACKSCLIPPERATTRIBUTES_H
#define CRACKSCLIPPERATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: CracksClipperAttributes
//
// Purpose:
//    Attributes for the cracks clipper operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class CracksClipperAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    CracksClipperAttributes();
    CracksClipperAttributes(const CracksClipperAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    CracksClipperAttributes(private_tmfs_t tmfs);
    CracksClipperAttributes(const CracksClipperAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~CracksClipperAttributes();

    virtual CracksClipperAttributes& operator = (const CracksClipperAttributes &obj);
    virtual bool operator == (const CracksClipperAttributes &obj) const;
    virtual bool operator != (const CracksClipperAttributes &obj) const;
private:
    void Init();
    void Copy(const CracksClipperAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectCrack1Var();
    void SelectCrack2Var();
    void SelectCrack3Var();
    void SelectStrainVar();
    void SelectInMassVar();

    // Property setting methods
    void SetCrack1Var(const std::string &crack1Var_);
    void SetCrack2Var(const std::string &crack2Var_);
    void SetCrack3Var(const std::string &crack3Var_);
    void SetStrainVar(const std::string &strainVar_);
    void SetShowCrack1(bool showCrack1_);
    void SetShowCrack2(bool showCrack2_);
    void SetShowCrack3(bool showCrack3_);
    void SetInMassVar(const std::string &inMassVar_);

    // Property getting methods
    const std::string &GetCrack1Var() const;
          std::string &GetCrack1Var();
    const std::string &GetCrack2Var() const;
          std::string &GetCrack2Var();
    const std::string &GetCrack3Var() const;
          std::string &GetCrack3Var();
    const std::string &GetStrainVar() const;
          std::string &GetStrainVar();
    bool              GetShowCrack1() const;
    bool              GetShowCrack2() const;
    bool              GetShowCrack3() const;
    const std::string &GetInMassVar() const;
          std::string &GetInMassVar();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const CracksClipperAttributes &) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_crack1Var = 0,
        ID_crack2Var,
        ID_crack3Var,
        ID_strainVar,
        ID_showCrack1,
        ID_showCrack2,
        ID_showCrack3,
        ID_inMassVar,
        ID__LAST
    };

private:
    std::string crack1Var;
    std::string crack2Var;
    std::string crack3Var;
    std::string strainVar;
    bool        showCrack1;
    bool        showCrack2;
    bool        showCrack3;
    std::string inMassVar;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define CRACKSCLIPPERATTRIBUTES_TMFS "ssssbbbs"

#endif
