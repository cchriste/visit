/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
// Creation:   omitted
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
        DomainNode,
        ZoneLabel,
        NodeLabel
    };
    enum CoordinateType
    {
        XY,
        RZ,
        ZR
    };
    enum TimeCurveType
    {
        Single_Y_Axis,
        Multiple_Y_Axes
    };

    // These constructors are for objects of this class
    PickAttributes();
    PickAttributes(const PickAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    PickAttributes(private_tmfs_t tmfs);
    PickAttributes(const PickAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~PickAttributes();

    virtual PickAttributes& operator = (const PickAttributes &obj);
    virtual bool operator == (const PickAttributes &obj) const;
    virtual bool operator != (const PickAttributes &obj) const;
private:
    void Init();
    void Copy(const PickAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVariables();
    void SelectPickLetter();
    void SelectIncidentElements();
    void SelectCellCoordinates();
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
    void SelectRangeOutput();
    void SelectElementLabel();
    void SelectSubsetName();
    void SelectFloatFormat();
    void SelectTimeOptions();
    void SelectPlotRequested();

    // Property setting methods
    void SetVariables(const stringVector &variables_);
    void SetShowIncidentElements(bool showIncidentElements_);
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
    void SetCellCoordinates(const doubleVector &cellCoordinates_);
    void SetTimeStep(int timeStep_);
    void SetDimension(int dimension_);
    void SetDatabaseName(const std::string &databaseName_);
    void SetActiveVariable(const std::string &activeVariable_);
    void SetPickPoint(const double *pickPoint_);
    void SetCellPoint(const double *cellPoint_);
    void SetNodePoint(const double *nodePoint_);
    void SetPlotBounds(const doubleVector &plotBounds_);
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
    void SetUseLabelAsPickLetter(bool useLabelAsPickLetter_);
    void SetShowGlobalIds(bool showGlobalIds_);
    void SetGlobalElement(int globalElement_);
    void SetGlobalIncidentElements(const intVector &globalIncidentElements_);
    void SetElementIsGlobal(bool elementIsGlobal_);
    void SetShowPickLetter(bool showPickLetter_);
    void SetHasRangeOutput(bool hasRangeOutput_);
    void SetRangeOutput(const MapNode &rangeOutput_);
    void SetElementLabel(const std::string &elementLabel_);
    void SetReusePickLetter(bool reusePickLetter_);
    void SetGhostType(int ghostType_);
    void SetHasMixedGhostTypes(int hasMixedGhostTypes_);
    void SetLinesData(bool linesData_);
    void SetShowPickHighlight(bool showPickHighlight_);
    void SetNotifyEnabled(bool notifyEnabled_);
    void SetInputTopoDim(int inputTopoDim_);
    void SetMeshCoordType(CoordinateType meshCoordType_);
    void SetCreateSpreadsheet(bool createSpreadsheet_);
    void SetSubsetName(const std::string &subsetName_);
    void SetFloatFormat(const std::string &floatFormat_);
    void SetTimePreserveCoord(bool timePreserveCoord_);
    void SetTimeCurveType(TimeCurveType timeCurveType_);
    void SetTimeOptions(const MapNode &timeOptions_);
    void SetPlotRequested(const MapNode &plotRequested_);

    // Property getting methods
    const stringVector &GetVariables() const;
          stringVector &GetVariables();
    bool               GetShowIncidentElements() const;
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
    const doubleVector &GetCellCoordinates() const;
          doubleVector &GetCellCoordinates();
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
    const doubleVector &GetPlotBounds() const;
          doubleVector &GetPlotBounds();
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
    bool               GetUseLabelAsPickLetter() const;
    bool               GetShowGlobalIds() const;
    int                GetGlobalElement() const;
    const intVector    &GetGlobalIncidentElements() const;
          intVector    &GetGlobalIncidentElements();
    bool               GetElementIsGlobal() const;
    bool               GetShowPickLetter() const;
    bool               GetHasRangeOutput() const;
    const MapNode      &GetRangeOutput() const;
          MapNode      &GetRangeOutput();
    const std::string  &GetElementLabel() const;
          std::string  &GetElementLabel();
    bool               GetReusePickLetter() const;
    int                GetGhostType() const;
    int                GetHasMixedGhostTypes() const;
    bool               GetLinesData() const;
    bool               GetShowPickHighlight() const;
    bool               GetNotifyEnabled() const;
    int                GetInputTopoDim() const;
    CoordinateType     GetMeshCoordType() const;
    bool               GetCreateSpreadsheet() const;
    const std::string  &GetSubsetName() const;
          std::string  &GetSubsetName();
    const std::string  &GetFloatFormat() const;
          std::string  &GetFloatFormat();
    bool               GetTimePreserveCoord() const;
    TimeCurveType      GetTimeCurveType() const;
    const MapNode      &GetTimeOptions() const;
          MapNode      &GetTimeOptions();
    const MapNode      &GetPlotRequested() const;
          MapNode      &GetPlotRequested();

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
    static std::string CoordinateType_ToString(CoordinateType);
    static bool CoordinateType_FromString(const std::string &, CoordinateType &);
protected:
    static std::string CoordinateType_ToString(int);
public:
    static std::string TimeCurveType_ToString(TimeCurveType);
    static bool TimeCurveType_FromString(const std::string &, TimeCurveType &);
protected:
    static std::string TimeCurveType_ToString(int);
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
    void AddLine(const double *_c0, const double *_c1, const int &pos);
    void Notify();
    void ClearLines();
    void CreateOutputMapNode(MapNode &m, bool withLetter);
    void CreateXMLString(std::string &os, bool withLetter = true);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_variables = 0,
        ID_showIncidentElements,
        ID_showNodeId,
        ID_showNodeDomainLogicalCoords,
        ID_showNodeBlockLogicalCoords,
        ID_showNodePhysicalCoords,
        ID_showZoneId,
        ID_showZoneDomainLogicalCoords,
        ID_showZoneBlockLogicalCoords,
        ID_clearWindow,
        ID_pickLetter,
        ID_fulfilled,
        ID_pickType,
        ID_domain,
        ID_elementNumber,
        ID_incidentElements,
        ID_cellCoordinates,
        ID_timeStep,
        ID_dimension,
        ID_databaseName,
        ID_activeVariable,
        ID_pickPoint,
        ID_cellPoint,
        ID_nodePoint,
        ID_plotBounds,
        ID_rayPoint1,
        ID_rayPoint2,
        ID_meshInfo,
        ID_realElementNumber,
        ID_realIncidentElements,
        ID_pnodeCoords,
        ID_dnodeCoords,
        ID_bnodeCoords,
        ID_dzoneCoords,
        ID_bzoneCoords,
        ID_needTransformMessage,
        ID_varInfo,
        ID_invalidVars,
        ID_doTimeCurve,
        ID_errorMessage,
        ID_error,
        ID_matSelected,
        ID_needActualCoords,
        ID_conciseOutput,
        ID_showTimeStep,
        ID_showMeshName,
        ID_blockPieceName,
        ID_groupPieceName,
        ID_ghosts,
        ID_includeGhosts,
        ID_elementIsGhost,
        ID_requiresGlyphPick,
        ID_locationSuccessful,
        ID_useLabelAsPickLetter,
        ID_showGlobalIds,
        ID_globalElement,
        ID_globalIncidentElements,
        ID_elementIsGlobal,
        ID_showPickLetter,
        ID_hasRangeOutput,
        ID_rangeOutput,
        ID_elementLabel,
        ID_reusePickLetter,
        ID_ghostType,
        ID_hasMixedGhostTypes,
        ID_linesData,
        ID_showPickHighlight,
        ID_notifyEnabled,
        ID_inputTopoDim,
        ID_meshCoordType,
        ID_createSpreadsheet,
        ID_subsetName,
        ID_floatFormat,
        ID_timePreserveCoord,
        ID_timeCurveType,
        ID_timeOptions,
        ID_plotRequested,
        ID__LAST
    };

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    stringVector         variables;
    bool                 showIncidentElements;
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
    doubleVector         cellCoordinates;
    int                  timeStep;
    int                  dimension;
    std::string          databaseName;
    std::string          activeVariable;
    double               pickPoint[3];
    double               cellPoint[3];
    double               nodePoint[3];
    doubleVector         plotBounds;
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
    bool                 useLabelAsPickLetter;
    bool                 showGlobalIds;
    int                  globalElement;
    intVector            globalIncidentElements;
    bool                 elementIsGlobal;
    bool                 showPickLetter;
    bool                 hasRangeOutput;
    MapNode              rangeOutput;
    std::string          elementLabel;
    bool                 reusePickLetter;
    int                  ghostType;
    int                  hasMixedGhostTypes;
    bool                 linesData;
    bool                 showPickHighlight;
    bool                 notifyEnabled;
    int                  inputTopoDim;
    int                  meshCoordType;
    bool                 createSpreadsheet;
    std::string          subsetName;
    std::string          floatFormat;
    bool                 timePreserveCoord;
    int                  timeCurveType;
    MapNode              timeOptions;
    MapNode              plotRequested;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define PICKATTRIBUTES_TMFS "s*bbbbbbbbbsbiiii*d*iissDDDd*DDsii*s*s*s*s*s*ba*s*bsbbbbbbssi*bbbbbbii*bbbmsbiibbbiibssbimm"

#endif
