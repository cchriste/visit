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

#ifndef FTLEATTRIBUTES_H
#define FTLEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: FTLEAttributes
//
// Purpose:
//    Attributes for FTLE
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class FTLEAttributes : public AttributeSubject
{
public:
    enum SourceType
    {
        NativeResolutionOfMesh,
        RegularGrid
    };
    enum Extents
    {
        Full,
        Subset
    };
    enum IntegrationDirection
    {
        Forward,
        Backward,
        Both
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
    enum SizeType
    {
        Absolute,
        FractionOfBBox
    };
    enum ParallelizationAlgorithmType
    {
        LoadOnDemand,
        ParallelStaticDomains,
        MasterSlave,
        VisItSelects
    };
    enum TerminationType
    {
        Time,
        Distance,
        Size
    };
    enum PathlinesCMFE
    {
        CONN_CMFE,
        POS_CMFE
    };

    // These constructors are for objects of this class
    FTLEAttributes();
    FTLEAttributes(const FTLEAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    FTLEAttributes(private_tmfs_t tmfs);
    FTLEAttributes(const FTLEAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~FTLEAttributes();

    virtual FTLEAttributes& operator = (const FTLEAttributes &obj);
    virtual bool operator == (const FTLEAttributes &obj) const;
    virtual bool operator != (const FTLEAttributes &obj) const;
private:
    void Init();
    void Copy(const FTLEAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectResolution();
    void SelectStartPosition();
    void SelectEndPosition();
    void SelectVelocitySource();

    // Property setting methods
    void SetSourceType(SourceType sourceType_);
    void SetResolution(const int *Resolution_);
    void SetUseDataSetStart(Extents UseDataSetStart_);
    void SetStartPosition(const double *StartPosition_);
    void SetUseDataSetEnd(Extents UseDataSetEnd_);
    void SetEndPosition(const double *EndPosition_);
    void SetIntegrationDirection(IntegrationDirection integrationDirection_);
    void SetMaxSteps(int maxSteps_);
    void SetTerminationType(TerminationType terminationType_);
    void SetTerminateBySize(bool terminateBySize_);
    void SetTermSize(double termSize_);
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
    void SetPathlinesCMFE(PathlinesCMFE pathlinesCMFE_);
    void SetForceNodeCenteredData(bool forceNodeCenteredData_);
    void SetIssueTerminationWarnings(bool issueTerminationWarnings_);
    void SetIssueStiffnessWarnings(bool issueStiffnessWarnings_);
    void SetIssueCriticalPointsWarnings(bool issueCriticalPointsWarnings_);
    void SetCriticalPointThreshold(double criticalPointThreshold_);

    // Property getting methods
    SourceType   GetSourceType() const;
    const int    *GetResolution() const;
          int    *GetResolution();
    Extents      GetUseDataSetStart() const;
    const double *GetStartPosition() const;
          double *GetStartPosition();
    Extents      GetUseDataSetEnd() const;
    const double *GetEndPosition() const;
          double *GetEndPosition();
    IntegrationDirection GetIntegrationDirection() const;
    int          GetMaxSteps() const;
    TerminationType GetTerminationType() const;
    bool         GetTerminateBySize() const;
    double       GetTermSize() const;
    bool         GetTerminateByDistance() const;
    double       GetTermDistance() const;
    bool         GetTerminateByTime() const;
    double       GetTermTime() const;
    double       GetMaxStepLength() const;
    bool         GetLimitMaximumTimestep() const;
    double       GetMaxTimeStep() const;
    double       GetRelTol() const;
    SizeType     GetAbsTolSizeType() const;
    double       GetAbsTolAbsolute() const;
    double       GetAbsTolBBox() const;
    FieldType    GetFieldType() const;
    double       GetFieldConstant() const;
    const double *GetVelocitySource() const;
          double *GetVelocitySource();
    IntegrationType GetIntegrationType() const;
    ParallelizationAlgorithmType GetParallelizationAlgorithmType() const;
    int          GetMaxProcessCount() const;
    int          GetMaxDomainCacheSize() const;
    int          GetWorkGroupSize() const;
    bool         GetPathlines() const;
    bool         GetPathlinesOverrideStartingTimeFlag() const;
    double       GetPathlinesOverrideStartingTime() const;
    PathlinesCMFE GetPathlinesCMFE() const;
    bool         GetForceNodeCenteredData() const;
    bool         GetIssueTerminationWarnings() const;
    bool         GetIssueStiffnessWarnings() const;
    bool         GetIssueCriticalPointsWarnings() const;
    double       GetCriticalPointThreshold() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string SourceType_ToString(SourceType);
    static bool SourceType_FromString(const std::string &, SourceType &);
protected:
    static std::string SourceType_ToString(int);
public:
    static std::string Extents_ToString(Extents);
    static bool Extents_FromString(const std::string &, Extents &);
protected:
    static std::string Extents_ToString(int);
public:
    static std::string IntegrationDirection_ToString(IntegrationDirection);
    static bool IntegrationDirection_FromString(const std::string &, IntegrationDirection &);
protected:
    static std::string IntegrationDirection_ToString(int);
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
    static std::string SizeType_ToString(SizeType);
    static bool SizeType_FromString(const std::string &, SizeType &);
protected:
    static std::string SizeType_ToString(int);
public:
    static std::string ParallelizationAlgorithmType_ToString(ParallelizationAlgorithmType);
    static bool ParallelizationAlgorithmType_FromString(const std::string &, ParallelizationAlgorithmType &);
protected:
    static std::string ParallelizationAlgorithmType_ToString(int);
public:
    static std::string TerminationType_ToString(TerminationType);
    static bool TerminationType_FromString(const std::string &, TerminationType &);
protected:
    static std::string TerminationType_ToString(int);
public:
    static std::string PathlinesCMFE_ToString(PathlinesCMFE);
    static bool PathlinesCMFE_FromString(const std::string &, PathlinesCMFE &);
protected:
    static std::string PathlinesCMFE_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const FTLEAttributes &) const;

    // IDs that can be used to identify fields in case statements
    enum {
        ID_sourceType = 0,
        ID_Resolution,
        ID_UseDataSetStart,
        ID_StartPosition,
        ID_UseDataSetEnd,
        ID_EndPosition,
        ID_integrationDirection,
        ID_maxSteps,
        ID_terminationType,
        ID_terminateBySize,
        ID_termSize,
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
        ID_pathlinesCMFE,
        ID_forceNodeCenteredData,
        ID_issueTerminationWarnings,
        ID_issueStiffnessWarnings,
        ID_issueCriticalPointsWarnings,
        ID_criticalPointThreshold,
        ID__LAST
    };

private:
    int    sourceType;
    int    Resolution[3];
    int    UseDataSetStart;
    double StartPosition[3];
    int    UseDataSetEnd;
    double EndPosition[3];
    int    integrationDirection;
    int    maxSteps;
    int    terminationType;
    bool   terminateBySize;
    double termSize;
    bool   terminateByDistance;
    double termDistance;
    bool   terminateByTime;
    double termTime;
    double maxStepLength;
    bool   limitMaximumTimestep;
    double maxTimeStep;
    double relTol;
    int    absTolSizeType;
    double absTolAbsolute;
    double absTolBBox;
    int    fieldType;
    double fieldConstant;
    double velocitySource[3];
    int    integrationType;
    int    parallelizationAlgorithmType;
    int    maxProcessCount;
    int    maxDomainCacheSize;
    int    workGroupSize;
    bool   pathlines;
    bool   pathlinesOverrideStartingTimeFlag;
    double pathlinesOverrideStartingTime;
    int    pathlinesCMFE;
    bool   forceNodeCenteredData;
    bool   issueTerminationWarnings;
    bool   issueStiffnessWarnings;
    bool   issueCriticalPointsWarnings;
    double criticalPointThreshold;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define FTLEATTRIBUTES_TMFS "iIiDiDiiibdbdbddbddiddidDiiiiibbdibbbbd"

#endif
