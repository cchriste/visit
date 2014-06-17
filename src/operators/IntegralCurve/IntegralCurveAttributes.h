/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

#ifndef INTEGRALCURVEATTRIBUTES_H
#define INTEGRALCURVEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: IntegralCurveAttributes
//
// Purpose:
//    Attributes for the IntegralCurve
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class IntegralCurveAttributes : public AttributeSubject
{
public:
    enum SourceType
    {
        Point,
        PointList,
        Line_,
        Circle,
        Plane,
        Sphere,
        Box,
        Selection
    };
    enum DataValue
    {
        Solid,
        SeedPointID,
        Speed,
        Vorticity,
        ArcLength,
        TimeAbsolute,
        TimeRelative,
        AverageDistanceFromSeed,
        CorrelationDistance,
        Difference,
        Variable
    };
    enum CropValue
    {
        Distance,
        Time,
        StepNumber
    };
    enum DisplayGeometry
    {
        Lines,
        Tubes,
        Ribbons
    };
    enum CoordinateSystem
    {
        AsIs,
        CylindricalToCartesian,
        CartesianToCylindrical
    };
    enum IntegrationDirection
    {
        Forward,
        Backward,
        Both
    };
    enum ParallelizationAlgorithmType
    {
        LoadOnDemand,
        ParallelStaticDomains,
        MasterSlave,
        VisItSelects
    };
    enum FieldType
    {
        Default,
        FlashField,
        M3DC12DField,
        M3DC13DField,
        Nek5000Field,
        NIMRODField
    };
    enum IntegrationType
    {
        Euler,
        Leapfrog,
        DormandPrince,
        AdamsBashforth,
        RK4,
        M3DC12DIntegrator
    };
    enum PathlinesCMFE
    {
        CONN_CMFE,
        POS_CMFE
    };
    enum SizeType
    {
        Absolute,
        FractionOfBBox
    };

    // These constructors are for objects of this class
    IntegralCurveAttributes();
    IntegralCurveAttributes(const IntegralCurveAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    IntegralCurveAttributes(private_tmfs_t tmfs);
    IntegralCurveAttributes(const IntegralCurveAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~IntegralCurveAttributes();

    virtual IntegralCurveAttributes& operator = (const IntegralCurveAttributes &obj);
    virtual bool operator == (const IntegralCurveAttributes &obj) const;
    virtual bool operator != (const IntegralCurveAttributes &obj) const;
private:
    void Init();
    void Copy(const IntegralCurveAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPointSource();
    void SelectLineStart();
    void SelectLineEnd();
    void SelectPlaneOrigin();
    void SelectPlaneNormal();
    void SelectPlaneUpAxis();
    void SelectSphereOrigin();
    void SelectBoxExtents();
    void SelectPointList();
    void SelectDataVariable();
    void SelectVelocitySource();
    void SelectSelection();

    // Property setting methods
    void SetSourceType(SourceType sourceType_);
    void SetPointSource(const double *pointSource_);
    void SetLineStart(const double *lineStart_);
    void SetLineEnd(const double *lineEnd_);
    void SetPlaneOrigin(const double *planeOrigin_);
    void SetPlaneNormal(const double *planeNormal_);
    void SetPlaneUpAxis(const double *planeUpAxis_);
    void SetRadius(double radius_);
    void SetSphereOrigin(const double *sphereOrigin_);
    void SetBoxExtents(const double *boxExtents_);
    void SetUseWholeBox(bool useWholeBox_);
    void SetPointList(const doubleVector &pointList_);
    void SetSampleDensity0(int sampleDensity0_);
    void SetSampleDensity1(int sampleDensity1_);
    void SetSampleDensity2(int sampleDensity2_);
    void SetDataValue(DataValue dataValue_);
    void SetDataVariable(const std::string &dataVariable_);
    void SetIntegrationDirection(IntegrationDirection integrationDirection_);
    void SetMaxSteps(int maxSteps_);
    void SetTerminateByDistance(bool terminateByDistance_);
    void SetTermDistance(double termDistance_);
    void SetTerminateByTime(bool terminateByTime_);
    void SetTermTime(double termTime_);
    void SetMaxStepLength(double maxStepLength_);
    void SetLimitMaximumTimestep(bool limitMaximumTimestep_);
    void SetMaxTimeStep(double maxTimeStep_);
    void SetRelTol(double relTol_);
    void SetAbsTolSizeType(SizeType absTolSizeType_);
    void SetAbsTolAbsolute(double absTolAbsolute_);
    void SetAbsTolBBox(double absTolBBox_);
    void SetFieldType(FieldType fieldType_);
    void SetFieldConstant(double fieldConstant_);
    void SetVelocitySource(const double *velocitySource_);
    void SetIntegrationType(IntegrationType integrationType_);
    void SetParallelizationAlgorithmType(ParallelizationAlgorithmType parallelizationAlgorithmType_);
    void SetMaxProcessCount(int maxProcessCount_);
    void SetMaxDomainCacheSize(int maxDomainCacheSize_);
    void SetWorkGroupSize(int workGroupSize_);
    void SetPathlines(bool pathlines_);
    void SetPathlinesOverrideStartingTimeFlag(bool pathlinesOverrideStartingTimeFlag_);
    void SetPathlinesOverrideStartingTime(double pathlinesOverrideStartingTime_);
    void SetPathlinesPeriod(double pathlinesPeriod_);
    void SetPathlinesCMFE(PathlinesCMFE pathlinesCMFE_);
    void SetDisplayGeometry(DisplayGeometry displayGeometry_);
    void SetCoordinateSystem(CoordinateSystem coordinateSystem_);
    void SetShowLines(bool showLines_);
    void SetShowPoints(bool showPoints_);
    void SetCropBeginFlag(bool cropBeginFlag_);
    void SetCropBegin(double cropBegin_);
    void SetCropEndFlag(bool cropEndFlag_);
    void SetCropEnd(double cropEnd_);
    void SetCropValue(CropValue cropValue_);
    void SetSampleDistance0(double sampleDistance0_);
    void SetSampleDistance1(double sampleDistance1_);
    void SetSampleDistance2(double sampleDistance2_);
    void SetFillInterior(bool fillInterior_);
    void SetRandomSamples(bool randomSamples_);
    void SetRandomSeed(int randomSeed_);
    void SetNumberOfRandomSamples(int numberOfRandomSamples_);
    void SetForceNodeCenteredData(bool forceNodeCenteredData_);
    void SetIssueTerminationWarnings(bool issueTerminationWarnings_);
    void SetIssueStiffnessWarnings(bool issueStiffnessWarnings_);
    void SetIssueCriticalPointsWarnings(bool issueCriticalPointsWarnings_);
    void SetCriticalPointThreshold(double criticalPointThreshold_);
    void SetCorrelationDistanceAngTol(double correlationDistanceAngTol_);
    void SetCorrelationDistanceMinDistAbsolute(double correlationDistanceMinDistAbsolute_);
    void SetCorrelationDistanceMinDistBBox(double correlationDistanceMinDistBBox_);
    void SetCorrelationDistanceMinDistType(SizeType correlationDistanceMinDistType_);
    void SetSelection(const std::string &selection_);

    // Property getting methods
    SourceType         GetSourceType() const;
    const double       *GetPointSource() const;
          double       *GetPointSource();
    const double       *GetLineStart() const;
          double       *GetLineStart();
    const double       *GetLineEnd() const;
          double       *GetLineEnd();
    const double       *GetPlaneOrigin() const;
          double       *GetPlaneOrigin();
    const double       *GetPlaneNormal() const;
          double       *GetPlaneNormal();
    const double       *GetPlaneUpAxis() const;
          double       *GetPlaneUpAxis();
    double             GetRadius() const;
    const double       *GetSphereOrigin() const;
          double       *GetSphereOrigin();
    const double       *GetBoxExtents() const;
          double       *GetBoxExtents();
    bool               GetUseWholeBox() const;
    const doubleVector &GetPointList() const;
          doubleVector &GetPointList();
    int                GetSampleDensity0() const;
    int                GetSampleDensity1() const;
    int                GetSampleDensity2() const;
    DataValue          GetDataValue() const;
    const std::string  &GetDataVariable() const;
          std::string  &GetDataVariable();
    IntegrationDirection GetIntegrationDirection() const;
    int                GetMaxSteps() const;
    bool               GetTerminateByDistance() const;
    double             GetTermDistance() const;
    bool               GetTerminateByTime() const;
    double             GetTermTime() const;
    double             GetMaxStepLength() const;
    bool               GetLimitMaximumTimestep() const;
    double             GetMaxTimeStep() const;
    double             GetRelTol() const;
    SizeType           GetAbsTolSizeType() const;
    double             GetAbsTolAbsolute() const;
    double             GetAbsTolBBox() const;
    FieldType          GetFieldType() const;
    double             GetFieldConstant() const;
    const double       *GetVelocitySource() const;
          double       *GetVelocitySource();
    IntegrationType    GetIntegrationType() const;
    ParallelizationAlgorithmType GetParallelizationAlgorithmType() const;
    int                GetMaxProcessCount() const;
    int                GetMaxDomainCacheSize() const;
    int                GetWorkGroupSize() const;
    bool               GetPathlines() const;
    bool               GetPathlinesOverrideStartingTimeFlag() const;
    double             GetPathlinesOverrideStartingTime() const;
    double             GetPathlinesPeriod() const;
    PathlinesCMFE      GetPathlinesCMFE() const;
    DisplayGeometry    GetDisplayGeometry() const;
    CoordinateSystem   GetCoordinateSystem() const;
    bool               GetShowLines() const;
    bool               GetShowPoints() const;
    bool               GetCropBeginFlag() const;
    double             GetCropBegin() const;
    bool               GetCropEndFlag() const;
    double             GetCropEnd() const;
    CropValue          GetCropValue() const;
    double             GetSampleDistance0() const;
    double             GetSampleDistance1() const;
    double             GetSampleDistance2() const;
    bool               GetFillInterior() const;
    bool               GetRandomSamples() const;
    int                GetRandomSeed() const;
    int                GetNumberOfRandomSamples() const;
    bool               GetForceNodeCenteredData() const;
    bool               GetIssueTerminationWarnings() const;
    bool               GetIssueStiffnessWarnings() const;
    bool               GetIssueCriticalPointsWarnings() const;
    double             GetCriticalPointThreshold() const;
    double             GetCorrelationDistanceAngTol() const;
    double             GetCorrelationDistanceMinDistAbsolute() const;
    double             GetCorrelationDistanceMinDistBBox() const;
    SizeType           GetCorrelationDistanceMinDistType() const;
    const std::string  &GetSelection() const;
          std::string  &GetSelection();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string SourceType_ToString(SourceType);
    static bool SourceType_FromString(const std::string &, SourceType &);
protected:
    static std::string SourceType_ToString(int);
public:
    static std::string DataValue_ToString(DataValue);
    static bool DataValue_FromString(const std::string &, DataValue &);
protected:
    static std::string DataValue_ToString(int);
public:
    static std::string CropValue_ToString(CropValue);
    static bool CropValue_FromString(const std::string &, CropValue &);
protected:
    static std::string CropValue_ToString(int);
public:
    static std::string DisplayGeometry_ToString(DisplayGeometry);
    static bool DisplayGeometry_FromString(const std::string &, DisplayGeometry &);
protected:
    static std::string DisplayGeometry_ToString(int);
public:
    static std::string CoordinateSystem_ToString(CoordinateSystem);
    static bool CoordinateSystem_FromString(const std::string &, CoordinateSystem &);
protected:
    static std::string CoordinateSystem_ToString(int);
public:
    static std::string IntegrationDirection_ToString(IntegrationDirection);
    static bool IntegrationDirection_FromString(const std::string &, IntegrationDirection &);
protected:
    static std::string IntegrationDirection_ToString(int);
public:
    static std::string ParallelizationAlgorithmType_ToString(ParallelizationAlgorithmType);
    static bool ParallelizationAlgorithmType_FromString(const std::string &, ParallelizationAlgorithmType &);
protected:
    static std::string ParallelizationAlgorithmType_ToString(int);
public:
    static std::string FieldType_ToString(FieldType);
    static bool FieldType_FromString(const std::string &, FieldType &);
protected:
    static std::string FieldType_ToString(int);
public:
    static std::string IntegrationType_ToString(IntegrationType);
    static bool IntegrationType_FromString(const std::string &, IntegrationType &);
protected:
    static std::string IntegrationType_ToString(int);
public:
    static std::string PathlinesCMFE_ToString(PathlinesCMFE);
    static bool PathlinesCMFE_FromString(const std::string &, PathlinesCMFE &);
protected:
    static std::string PathlinesCMFE_ToString(int);
public:
    static std::string SizeType_ToString(SizeType);
    static bool SizeType_FromString(const std::string &, SizeType &);
protected:
    static std::string SizeType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const IntegralCurveAttributes &) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_sourceType = 0,
        ID_pointSource,
        ID_lineStart,
        ID_lineEnd,
        ID_planeOrigin,
        ID_planeNormal,
        ID_planeUpAxis,
        ID_radius,
        ID_sphereOrigin,
        ID_boxExtents,
        ID_useWholeBox,
        ID_pointList,
        ID_sampleDensity0,
        ID_sampleDensity1,
        ID_sampleDensity2,
        ID_dataValue,
        ID_dataVariable,
        ID_integrationDirection,
        ID_maxSteps,
        ID_terminateByDistance,
        ID_termDistance,
        ID_terminateByTime,
        ID_termTime,
        ID_maxStepLength,
        ID_limitMaximumTimestep,
        ID_maxTimeStep,
        ID_relTol,
        ID_absTolSizeType,
        ID_absTolAbsolute,
        ID_absTolBBox,
        ID_fieldType,
        ID_fieldConstant,
        ID_velocitySource,
        ID_integrationType,
        ID_parallelizationAlgorithmType,
        ID_maxProcessCount,
        ID_maxDomainCacheSize,
        ID_workGroupSize,
        ID_pathlines,
        ID_pathlinesOverrideStartingTimeFlag,
        ID_pathlinesOverrideStartingTime,
        ID_pathlinesPeriod,
        ID_pathlinesCMFE,
        ID_displayGeometry,
        ID_coordinateSystem,
        ID_showLines,
        ID_showPoints,
        ID_cropBeginFlag,
        ID_cropBegin,
        ID_cropEndFlag,
        ID_cropEnd,
        ID_cropValue,
        ID_sampleDistance0,
        ID_sampleDistance1,
        ID_sampleDistance2,
        ID_fillInterior,
        ID_randomSamples,
        ID_randomSeed,
        ID_numberOfRandomSamples,
        ID_forceNodeCenteredData,
        ID_issueTerminationWarnings,
        ID_issueStiffnessWarnings,
        ID_issueCriticalPointsWarnings,
        ID_criticalPointThreshold,
        ID_correlationDistanceAngTol,
        ID_correlationDistanceMinDistAbsolute,
        ID_correlationDistanceMinDistBBox,
        ID_correlationDistanceMinDistType,
        ID_selection,
        ID__LAST
    };

private:
    int          sourceType;
    double       pointSource[3];
    double       lineStart[3];
    double       lineEnd[3];
    double       planeOrigin[3];
    double       planeNormal[3];
    double       planeUpAxis[3];
    double       radius;
    double       sphereOrigin[3];
    double       boxExtents[6];
    bool         useWholeBox;
    doubleVector pointList;
    int          sampleDensity0;
    int          sampleDensity1;
    int          sampleDensity2;
    int          dataValue;
    std::string  dataVariable;
    int          integrationDirection;
    int          maxSteps;
    bool         terminateByDistance;
    double       termDistance;
    bool         terminateByTime;
    double       termTime;
    double       maxStepLength;
    bool         limitMaximumTimestep;
    double       maxTimeStep;
    double       relTol;
    int          absTolSizeType;
    double       absTolAbsolute;
    double       absTolBBox;
    int          fieldType;
    double       fieldConstant;
    double       velocitySource[3];
    int          integrationType;
    int          parallelizationAlgorithmType;
    int          maxProcessCount;
    int          maxDomainCacheSize;
    int          workGroupSize;
    bool         pathlines;
    bool         pathlinesOverrideStartingTimeFlag;
    double       pathlinesOverrideStartingTime;
    double       pathlinesPeriod;
    int          pathlinesCMFE;
    int          displayGeometry;
    int          coordinateSystem;
    bool         showLines;
    bool         showPoints;
    bool         cropBeginFlag;
    double       cropBegin;
    bool         cropEndFlag;
    double       cropEnd;
    int          cropValue;
    double       sampleDistance0;
    double       sampleDistance1;
    double       sampleDistance2;
    bool         fillInterior;
    bool         randomSamples;
    int          randomSeed;
    int          numberOfRandomSamples;
    bool         forceNodeCenteredData;
    bool         issueTerminationWarnings;
    bool         issueStiffnessWarnings;
    bool         issueCriticalPointsWarnings;
    double       criticalPointThreshold;
    double       correlationDistanceAngTol;
    double       correlationDistanceMinDistAbsolute;
    double       correlationDistanceMinDistBBox;
    int          correlationDistanceMinDistType;
    std::string  selection;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define INTEGRALCURVEATTRIBUTES_TMFS "iDDDDDDdDDbd*iiiisiibdbddbddiddidDiiiiibbddiiibbbdbdidddbbiibbbbddddis"

#endif
