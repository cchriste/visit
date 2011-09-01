/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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

#ifndef QUERYLIST_H
#define QUERYLIST_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: QueryList
//
// Purpose:
//    List of supported queries
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API QueryList : public AttributeSubject
{
public:
    enum QueryType
    {
        DatabaseQuery,
        PointQuery,
        LineQuery
    };
    enum WindowType
    {
        Basic,
        DomainNode,
        DomainNodeVars,
        DomainZone,
        DomainZoneVars,
        ActualData,
        ActualDataVars,
        LineDistribution,
        HohlraumFlux,
        ConnCompSummary,
        ShapeletsDecomp,
        XRayImage,
        StreamlineInfo,
        Pick,
        Lineout
    };
    enum Groups
    {
        CurveRelated,
        MeshRelated,
        TimeRelated,
        VariableRelated,
        ShapeRelated,
        ConnectedComponentsRelated,
        Miscellaneous,
        NumGroups
    };
    enum QueryMode
    {
        QueryOnly,
        QueryAndTime,
        TimeOnly
    };

    // These constructors are for objects of this class
    QueryList();
    QueryList(const QueryList &obj);
protected:
    // These constructors are for objects derived from this class
    QueryList(private_tmfs_t tmfs);
    QueryList(const QueryList &obj, private_tmfs_t tmfs);
public:
    virtual ~QueryList();

    virtual QueryList& operator = (const QueryList &obj);
    virtual bool operator == (const QueryList &obj) const;
    virtual bool operator != (const QueryList &obj) const;
private:
    void Init();
    void Copy(const QueryList &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectNames();
    void SelectTypes();
    void SelectGroups();
    void SelectNumInputs();
    void SelectAllowedVarTypes();
    void SelectWinType();
    void SelectQueryMode();
    void SelectNumVars();
    void SelectCanBePublic();
    void SelectRequiresVarSelection();

    // Property setting methods
    void SetNames(const stringVector &names_);
    void SetTypes(const intVector &types_);
    void SetGroups(const intVector &groups_);
    void SetNumInputs(const intVector &numInputs_);
    void SetAllowedVarTypes(const intVector &allowedVarTypes_);
    void SetWinType(const intVector &winType_);
    void SetQueryMode(const intVector &queryMode_);
    void SetNumVars(const intVector &numVars_);
    void SetCanBePublic(const intVector &canBePublic_);
    void SetRequiresVarSelection(const intVector &requiresVarSelection_);

    // Property getting methods
    const stringVector &GetNames() const;
          stringVector &GetNames();
    const intVector    &GetTypes() const;
          intVector    &GetTypes();
    const intVector    &GetGroups() const;
          intVector    &GetGroups();
    const intVector    &GetNumInputs() const;
          intVector    &GetNumInputs();
    const intVector    &GetAllowedVarTypes() const;
          intVector    &GetAllowedVarTypes();
    const intVector    &GetWinType() const;
          intVector    &GetWinType();
    const intVector    &GetQueryMode() const;
          intVector    &GetQueryMode();
    const intVector    &GetNumVars() const;
          intVector    &GetNumVars();
    const intVector    &GetCanBePublic() const;
          intVector    &GetCanBePublic();
    const intVector    &GetRequiresVarSelection() const;
          intVector    &GetRequiresVarSelection();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string QueryType_ToString(QueryType);
    static bool QueryType_FromString(const std::string &, QueryType &);
protected:
    static std::string QueryType_ToString(int);
public:
    static std::string WindowType_ToString(WindowType);
    static bool WindowType_FromString(const std::string &, WindowType &);
protected:
    static std::string WindowType_ToString(int);
public:
    static std::string Groups_ToString(Groups);
    static bool Groups_FromString(const std::string &, Groups &);
protected:
    static std::string Groups_ToString(int);
public:
    static std::string QueryMode_ToString(QueryMode);
    static bool QueryMode_FromString(const std::string &, QueryMode &);
protected:
    static std::string QueryMode_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void AddQuery(const std::string &name, QueryType t, Groups g, WindowType w, int num_input, int allowedVars, QueryMode qMode, int numVars = 1, int reqVars = 0);
    bool QueryExists(const std::string &name, QueryType t);
    int NumberOfInputsForQuery(const std::string &name);
    int AllowedVarsForQuery(const std::string &name);
    bool TimeQueryAvailable(const std::string &name) ;
    int GetWindowType(const std::string &name) ;
    int NumberOfVarsForQuery(const std::string &name);
    bool RegularQueryAvailable(const std::string &name) ;
    int GetQueryType(const std::string &name) ;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_names = 0,
        ID_types,
        ID_groups,
        ID_numInputs,
        ID_allowedVarTypes,
        ID_winType,
        ID_queryMode,
        ID_numVars,
        ID_canBePublic,
        ID_requiresVarSelection,
        ID__LAST
    };

private:
    stringVector names;
    intVector    types;
    intVector    groups;
    intVector    numInputs;
    intVector    allowedVarTypes;
    intVector    winType;
    intVector    queryMode;
    intVector    numVars;
    intVector    canBePublic;
    intVector    requiresVarSelection;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define QUERYLIST_TMFS "s*i*i*i*i*i*i*i*i*i*"

#endif
