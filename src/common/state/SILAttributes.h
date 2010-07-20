/*****************************************************************************
*
* Copyright (c) 2000 - 2010, Lawrence Livermore National Security, LLC
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

#ifndef SILATTRIBUTES_H
#define SILATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
class NamespaceAttributes;
class SILMatrixAttributes;
class SILArrayAttributes;

// ****************************************************************************
// Class: SILAttributes
//
// Purpose:
//    This class contains the information needed to represent a SIL.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SILAttributes : public AttributeSubject
{
public:
    SILAttributes();
    SILAttributes(const SILAttributes &obj);
    virtual ~SILAttributes();

    virtual SILAttributes& operator = (const SILAttributes &obj);
    virtual bool operator == (const SILAttributes &obj) const;
    virtual bool operator != (const SILAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSetNames();
    void SelectSetIds();
    void SelectWholeList();
    void SelectCategory();
    void SelectRole();
    void SelectSuperset();
    void SelectNspace();
    void SelectMatrices();
    void SelectArrays();
    void SelectOrder();

    // Property setting methods
    void SetNSets(int nSets_);
    void SetSetNames(const stringVector &setNames_);
    void SetSetIds(const intVector &setIds_);
    void SetWholeList(const intVector &wholeList_);
    void SetNCollections(int nCollections_);
    void SetCategory(const stringVector &category_);
    void SetRole(const intVector &role_);
    void SetSuperset(const intVector &superset_);
    void SetOrder(const intVector &order_);

    // Property getting methods
    int                GetNSets() const;
    const stringVector &GetSetNames() const;
          stringVector &GetSetNames();
    const intVector    &GetSetIds() const;
          intVector    &GetSetIds();
    const intVector    &GetWholeList() const;
          intVector    &GetWholeList();
    int                GetNCollections() const;
    const stringVector &GetCategory() const;
          stringVector &GetCategory();
    const intVector    &GetRole() const;
          intVector    &GetRole();
    const intVector    &GetSuperset() const;
          intVector    &GetSuperset();
    const AttributeGroupVector &GetNspace() const;
          AttributeGroupVector &GetNspace();
    const AttributeGroupVector &GetMatrices() const;
          AttributeGroupVector &GetMatrices();
    const AttributeGroupVector &GetArrays() const;
          AttributeGroupVector &GetArrays();
    const intVector    &GetOrder() const;
          intVector    &GetOrder();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddNspace(const NamespaceAttributes &);
    void ClearNspaces();
    void RemoveNspace(int i);
    int  GetNumNspaces() const;
    NamespaceAttributes &GetNspace(int i);
    const NamespaceAttributes &GetNspace(int i) const;

    void AddMatrices(const SILMatrixAttributes &);
    void ClearMatrices();
    void RemoveMatrices(int i);
    int  GetNumMatrices() const;
    SILMatrixAttributes &GetMatrices(int i);
    const SILMatrixAttributes &GetMatrices(int i) const;

    void AddArrays(const SILArrayAttributes &);
    void ClearArrays();
    void RemoveArrays(int i);
    int  GetNumArrays() const;
    SILArrayAttributes &GetArrays(int i);
    const SILArrayAttributes &GetArrays(int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_nSets = 0,
        ID_setNames,
        ID_setIds,
        ID_wholeList,
        ID_nCollections,
        ID_category,
        ID_role,
        ID_superset,
        ID_nspace,
        ID_matrices,
        ID_arrays,
        ID_order
    };

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    int                  nSets;
    stringVector         setNames;
    intVector            setIds;
    intVector            wholeList;
    int                  nCollections;
    stringVector         category;
    intVector            role;
    intVector            superset;
    AttributeGroupVector nspace;
    AttributeGroupVector matrices;
    AttributeGroupVector arrays;
    intVector            order;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
