/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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

#ifndef PICKATTRIBUTES_H
#define PICKATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
class PickVarInfo;
#include <visitstream.h>

// ****************************************************************************
// Class: PickAttributes
//
// Purpose:
//    This class contains attributes used for pick.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Nov 28 14:13:09 PST 2007
//
// Modifications:
//   
// ****************************************************************************

class STATE_API PickAttributes : public AttributeSubject
{
public:
    enum PickType
    {
        Zone,
        Node,
        CurveZone,
        CurveNode,
        DomainZone,
        DomainNode
    };

    PickAttributes();
    PickAttributes(const PickAttributes &obj);
    virtual ~PickAttributes();

    virtual PickAttributes& operator = (const PickAttributes &obj);
    virtual bool operator == (const PickAttributes &obj) const;
    virtual bool operator != (const PickAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVariables();
    void SelectPickLetter();
    void SelectIncidentElements();
    void SelectDatabaseName();
    void SelectActiveVariable();
    void SelectPickPoint();
    void SelectCellPoint();
    void SelectNodePoint();
    void SelectPlotBounds();
    void SelectRayPoint1();
    void SelectRayPoint2();
    void SelectMeshInfo();
    void SelectRealIncidentElements();
    void SelectPnodeCoords();
    void SelectDnodeCoords();
    void SelectBnodeCoords();
    void SelectDzoneCoords();
    void SelectBzoneCoords();
    void SelectVarInfo();
    void SelectInvalidVars();
    void SelectErrorMessage();
    void SelectBlockPieceName();
    void SelectGroupPieceName();
    void SelectGhosts();
    void SelectGlobalIncidentElements();
    void SelectSubsetName();
    void SelectFloatFormat();

    // Property setting methods
    void SetVariables(const stringVector &variables_);
    void SetDisplayIncidentElements(bool displayIncidentElements_);
    void SetShowNodeId(bool showNodeId_);
    void SetShowNodeDomainLogicalCoords(bool showNodeDomainLogicalCoords_);
    void SetShowNodeBlockLogicalCoords(bool showNodeBlockLogicalCoords_);
    void SetShowNodePhysicalCoords(bool showNodePhysicalCoords_);
    void SetShowZoneId(bool showZoneId_);
    void SetShowZoneDomainLogicalCoords(bool showZoneDomainLogicalCoords_);
    void SetShowZoneBlockLogicalCoords(bool showZoneBlockLogicalCoords_);
    void SetClearWindow(bool clearWindow_);
    void SetPickLetter(const std::string &pickLetter_);
    void SetFulfilled(bool fulfilled_);
    void SetPickType(PickType pickType_);
    void SetDomain(int domain_);
    void SetElementNumber(int elementNumber_);
    void SetIncidentElements(const intVector &incidentElements_);
    void SetTimeStep(int timeStep_);
    void SetDimension(int dimension_);
    void SetDatabaseName(const std::string &databaseName_);
    void SetActiveVariable(const std::string &activeVariable_);
    void SetPickPoint(const double *pickPoint_);
    void SetCellPoint(const double *cellPoint_);
    void SetNodePoint(const double *nodePoint_);
    void SetPlotBounds(const double *plotBounds_);
    void SetRayPoint1(const double *rayPoint1_);
    void SetRayPoint2(const double *rayPoint2_);
    void SetMeshInfo(const std::string &meshInfo_);
    void SetRealElementNumber(int realElementNumber_);
    void SetRealIncidentElements(const intVector &realIncidentElements_);
    void SetPnodeCoords(const stringVector &pnodeCoords_);
    void SetDnodeCoords(const stringVector &dnodeCoords_);
    void SetBnodeCoords(const stringVector &bnodeCoords_);
    void SetDzoneCoords(const stringVector &dzoneCoords_);
    void SetBzoneCoords(const stringVector &bzoneCoords_);
    void SetNeedTransformMessage(bool needTransformMessage_);
    void SetInvalidVars(const stringVector &invalidVars_);
    void SetDoTimeCurve(bool doTimeCurve_);
    void SetErrorMessage(const std::string &errorMessage_);
    void SetError(bool error_);
    void SetMatSelected(bool matSelected_);
    void SetNeedActualCoords(bool needActualCoords_);
    void SetConciseOutput(bool conciseOutput_);
    void SetShowTimeStep(bool showTimeStep_);
    void SetShowMeshName(bool showMeshName_);
    void SetBlockPieceName(const std::string &blockPieceName_);
    void SetGroupPieceName(const std::string &groupPieceName_);
    void SetGhosts(const intVector &ghosts_);
    void SetIncludeGhosts(bool includeGhosts_);
    void SetElementIsGhost(bool elementIsGhost_);
    void SetRequiresGlyphPick(bool requiresGlyphPick_);
    void SetLocationSuccessful(bool locationSuccessful_);
    void SetDisplayGlobalIds(bool displayGlobalIds_);
    void SetGlobalElement(int globalElement_);
    void SetGlobalIncidentElements(const intVector &globalIncidentElements_);
    void SetElementIsGlobal(bool elementIsGlobal_);
    void SetDisplayPickLetter(bool displayPickLetter_);
    void SetGhostType(int ghostType_);
    void SetHasMixedGhostTypes(int hasMixedGhostTypes_);
    void SetLinesData(bool linesData_);
    void SetInputTopoDim(int inputTopoDim_);
    void SetMeshCoordType(int meshCoordType_);
    void SetCreateSpreadsheet(bool createSpreadsheet_);
    void SetSubsetName(const std::string &subsetName_);
    void SetFloatFormat(const std::string &floatFormat_);
    void SetTimePreserveCoord(bool timePreserveCoord_);

    // Property getting methods
    const stringVector &GetVariables() const;
          stringVector &GetVariables();
    bool               GetDisplayIncidentElements() const;
    bool               GetShowNodeId() const;
    bool               GetShowNodeDomainLogicalCoords() const;
    bool               GetShowNodeBlockLogicalCoords() const;
    bool               GetShowNodePhysicalCoords() const;
    bool               GetShowZoneId() const;
    bool               GetShowZoneDomainLogicalCoords() const;
    bool               GetShowZoneBlockLogicalCoords() const;
    bool               GetClearWindow() const;
    const std::string  &GetPickLetter() const;
          std::string  &GetPickLetter();
    bool               GetFulfilled() const;
    PickType           GetPickType() const;
    int                GetDomain() const;
    int                GetElementNumber() const;
    const intVector    &GetIncidentElements() const;
          intVector    &GetIncidentElements();
    int                GetTimeStep() const;
    int                GetDimension() const;
    const std::string  &GetDatabaseName() const;
          std::string  &GetDatabaseName();
    const std::string  &GetActiveVariable() const;
          std::string  &GetActiveVariable();
    const double       *GetPickPoint() const;
          double       *GetPickPoint();
    const double       *GetCellPoint() const;
          double       *GetCellPoint();
    const double       *GetNodePoint() const;
          double       *GetNodePoint();
    const double       *GetPlotBounds() const;
          double       *GetPlotBounds();
    const double       *GetRayPoint1() const;
          double       *GetRayPoint1();
    const double       *GetRayPoint2() const;
          double       *GetRayPoint2();
    const std::string  &GetMeshInfo() const;
          std::string  &GetMeshInfo();
    int                GetRealElementNumber() const;
    const intVector    &GetRealIncidentElements() const;
          intVector    &GetRealIncidentElements();
    const stringVector &GetPnodeCoords() const;
          stringVector &GetPnodeCoords();
    const stringVector &GetDnodeCoords() const;
          stringVector &GetDnodeCoords();
    const stringVector &GetBnodeCoords() const;
          stringVector &GetBnodeCoords();
    const stringVector &GetDzoneCoords() const;
          stringVector &GetDzoneCoords();
    const stringVector &GetBzoneCoords() const;
          stringVector &GetBzoneCoords();
    bool               GetNeedTransformMessage() const;
    const AttributeGroupVector &GetVarInfo() const;
          AttributeGroupVector &GetVarInfo();
    const stringVector &GetInvalidVars() const;
          stringVector &GetInvalidVars();
    bool               GetDoTimeCurve() const;
    const std::string  &GetErrorMessage() const;
          std::string  &GetErrorMessage();
    bool               GetError() const;
    bool               GetMatSelected() const;
    bool               GetNeedActualCoords() const;
    bool               GetConciseOutput() const;
    bool               GetShowTimeStep() const;
    bool               GetShowMeshName() const;
    const std::string  &GetBlockPieceName() const;
          std::string  &GetBlockPieceName();
    const std::string  &GetGroupPieceName() const;
          std::string  &GetGroupPieceName();
    const intVector    &GetGhosts() const;
          intVector    &GetGhosts();
    bool               GetIncludeGhosts() const;
    bool               GetElementIsGhost() const;
    bool               GetRequiresGlyphPick() const;
    bool               GetLocationSuccessful() const;
    bool               GetDisplayGlobalIds() const;
    int                GetGlobalElement() const;
    const intVector    &GetGlobalIncidentElements() const;
          intVector    &GetGlobalIncidentElements();
    bool               GetElementIsGlobal() const;
    bool               GetDisplayPickLetter() const;
    int                GetGhostType() const;
    int                GetHasMixedGhostTypes() const;
    bool               GetLinesData() const;
    int                GetInputTopoDim() const;
    int                GetMeshCoordType() const;
    bool               GetCreateSpreadsheet() const;
    const std::string  &GetSubsetName() const;
          std::string  &GetSubsetName();
    const std::string  &GetFloatFormat() const;
          std::string  &GetFloatFormat();
    bool               GetTimePreserveCoord() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddVarInfo(const PickVarInfo &);
    void ClearVarInfos();
    void RemoveVarInfo(int i);
    int  GetNumVarInfos() const;
    PickVarInfo &GetVarInfo(int i);
    const PickVarInfo &GetVarInfo(int i) const;

    PickVarInfo &operator [] (int i);
    const PickVarInfo &operator [] (int i) const;

    // Enum conversion functions
    static std::string PickType_ToString(PickType);
    static bool PickType_FromString(const std::string &, PickType &);
protected:
    static std::string PickType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void PrintSelf(ostream &os);
    void CreateOutputString(std::string &os, bool withLetter = true);
    void PrepareForNewPick();
    void CreateConciseOutputString(std::string &os, bool withLetter = true);
    void SetRayPoint1(const doubleVector &);
    void SetRayPoint2(const doubleVector &);
protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    stringVector         variables;
    bool                 displayIncidentElements;
    bool                 showNodeId;
    bool                 showNodeDomainLogicalCoords;
    bool                 showNodeBlockLogicalCoords;
    bool                 showNodePhysicalCoords;
    bool                 showZoneId;
    bool                 showZoneDomainLogicalCoords;
    bool                 showZoneBlockLogicalCoords;
    bool                 clearWindow;
    std::string          pickLetter;
    bool                 fulfilled;
    int                  pickType;
    int                  domain;
    int                  elementNumber;
    intVector            incidentElements;
    int                  timeStep;
    int                  dimension;
    std::string          databaseName;
    std::string          activeVariable;
    double               pickPoint[3];
    double               cellPoint[3];
    double               nodePoint[3];
    double               plotBounds[6];
    double               rayPoint1[3];
    double               rayPoint2[3];
    std::string          meshInfo;
    int                  realElementNumber;
    intVector            realIncidentElements;
    stringVector         pnodeCoords;
    stringVector         dnodeCoords;
    stringVector         bnodeCoords;
    stringVector         dzoneCoords;
    stringVector         bzoneCoords;
    bool                 needTransformMessage;
    AttributeGroupVector varInfo;
    stringVector         invalidVars;
    bool                 doTimeCurve;
    std::string          errorMessage;
    bool                 error;
    bool                 matSelected;
    bool                 needActualCoords;
    bool                 conciseOutput;
    bool                 showTimeStep;
    bool                 showMeshName;
    std::string          blockPieceName;
    std::string          groupPieceName;
    intVector            ghosts;
    bool                 includeGhosts;
    bool                 elementIsGhost;
    bool                 requiresGlyphPick;
    bool                 locationSuccessful;
    bool                 displayGlobalIds;
    int                  globalElement;
    intVector            globalIncidentElements;
    bool                 elementIsGlobal;
    bool                 displayPickLetter;
    int                  ghostType;
    int                  hasMixedGhostTypes;
    bool                 linesData;
    int                  inputTopoDim;
    int                  meshCoordType;
    bool                 createSpreadsheet;
    std::string          subsetName;
    std::string          floatFormat;
    bool                 timePreserveCoord;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
