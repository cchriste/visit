/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

#ifndef FIVEFOLDTETSUBDIVISIONATTRIBUTES_H
#define FIVEFOLDTETSUBDIVISIONATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: FiveFoldTetSubdivisionAttributes
//
// Purpose:
//    Attributes for five fold tetrahedral subdivision operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class FiveFoldTetSubdivisionAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    FiveFoldTetSubdivisionAttributes();
    FiveFoldTetSubdivisionAttributes(const FiveFoldTetSubdivisionAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    FiveFoldTetSubdivisionAttributes(private_tmfs_t tmfs);
    FiveFoldTetSubdivisionAttributes(const FiveFoldTetSubdivisionAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~FiveFoldTetSubdivisionAttributes();

    virtual FiveFoldTetSubdivisionAttributes& operator = (const FiveFoldTetSubdivisionAttributes &obj);
    virtual bool operator == (const FiveFoldTetSubdivisionAttributes &obj) const;
    virtual bool operator != (const FiveFoldTetSubdivisionAttributes &obj) const;
private:
    void Init();
    void Copy(const FiveFoldTetSubdivisionAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectIdVar();
    void SelectValueVar();
    void SelectContourTreeFilename();
    void SelectSelectedIds();
    void SelectHighlightedIds();

    // Property setting methods
    void SetOddParityHasSixNeighborhood(bool oddParityHasSixNeighborhood_);
    void SetAddComponentInformation(bool addComponentInformation_);
    void SetIdVar(const std::string &idVar_);
    void SetValueVar(const std::string &valueVar_);
    void SetContourTreeFilename(const std::string &contourTreeFilename_);
    void SetIsovalue(double isovalue_);
    void SetSelectedIds(const intVector &selectedIds_);
    void SetHighlightedIds(const intVector &highlightedIds_);

    // Property getting methods
    bool              GetOddParityHasSixNeighborhood() const;
    bool              GetAddComponentInformation() const;
    const std::string &GetIdVar() const;
          std::string &GetIdVar();
    const std::string &GetValueVar() const;
          std::string &GetValueVar();
    const std::string &GetContourTreeFilename() const;
          std::string &GetContourTreeFilename();
    double            GetIsovalue() const;
    const intVector   &GetSelectedIds() const;
          intVector   &GetSelectedIds();
    const intVector   &GetHighlightedIds() const;
          intVector   &GetHighlightedIds();

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
        ID_oddParityHasSixNeighborhood = 0,
        ID_addComponentInformation,
        ID_idVar,
        ID_valueVar,
        ID_contourTreeFilename,
        ID_isovalue,
        ID_selectedIds,
        ID_highlightedIds,
        ID__LAST
    };

private:
    bool        oddParityHasSixNeighborhood;
    bool        addComponentInformation;
    std::string idVar;
    std::string valueVar;
    std::string contourTreeFilename;
    double      isovalue;
    intVector   selectedIds;
    intVector   highlightedIds;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define FIVEFOLDTETSUBDIVISIONATTRIBUTES_TMFS "bbsssdi*i*"

#endif
