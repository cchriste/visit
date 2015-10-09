/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

#ifndef SELECTIONLIST_H
#define SELECTIONLIST_H
#include <state_exports.h>
#include <AttributeSubject.h>

class SelectionProperties;
class SelectionSummary;

// ****************************************************************************
// Class: SelectionList
//
// Purpose:
//    Contains a list of SelectionProperties objects.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SelectionList : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    SelectionList();
    SelectionList(const SelectionList &obj);
protected:
    // These constructors are for objects derived from this class
    SelectionList(private_tmfs_t tmfs);
    SelectionList(const SelectionList &obj, private_tmfs_t tmfs);
public:
    virtual ~SelectionList();

    virtual SelectionList& operator = (const SelectionList &obj);
    virtual bool operator == (const SelectionList &obj) const;
    virtual bool operator != (const SelectionList &obj) const;
private:
    void Init();
    void Copy(const SelectionList &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSelections();
    void SelectSelectionSummary();

    // Property setting methods
    void SetAutoApplyUpdates(bool autoApplyUpdates_);

    // Property getting methods
    bool GetAutoApplyUpdates() const;
    const AttributeGroupVector &GetSelections() const;
          AttributeGroupVector &GetSelections();
    const AttributeGroupVector &GetSelectionSummary() const;
          AttributeGroupVector &GetSelectionSummary();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddSelections(const SelectionProperties &);
    void ClearSelections();
    void RemoveSelections(int i);
    int  GetNumSelections() const;
    SelectionProperties &GetSelections(int i);
    const SelectionProperties &GetSelections(int i) const;

    void AddSelectionSummary(const SelectionSummary &);
    void ClearSelectionSummarys();
    void RemoveSelectionSummary(int i);
    int  GetNumSelectionSummarys() const;
    SelectionSummary &GetSelectionSummary(int i);
    const SelectionSummary &GetSelectionSummary(int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    int GetSelection(const std::string &name) const;
    int GetSelectionSummary(const std::string &name) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_autoApplyUpdates = 0,
        ID_selections,
        ID_selectionSummary,
        ID__LAST
    };

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    bool                 autoApplyUpdates;
    AttributeGroupVector selections;
    AttributeGroupVector selectionSummary;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define SELECTIONLIST_TMFS "ba*a*"

#endif
