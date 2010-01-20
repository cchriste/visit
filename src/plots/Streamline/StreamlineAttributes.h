/*****************************************************************************
*
* Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
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

#ifndef STREAMLINEATTRIBUTES_H
#define STREAMLINEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>

// ****************************************************************************
// Class: StreamlineAttributes
//
// Purpose:
//    Attributes for the Streamline plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class StreamlineAttributes : public AttributeSubject
{
public:
    enum SourceType
    {
        SpecifiedPoint,
        SpecifiedLine,
        SpecifiedPlane,
        SpecifiedSphere,
        SpecifiedBox,
        SpecifiedCircle,
        SpecifiedPointList
    };
    enum ColoringMethod
    {
        Solid,
        ColorBySpeed,
        ColorByVorticity,
        ColorByLength,
        ColorByTime,
        ColorBySeedPointID,
        ColorByVariable
    };
    enum DisplayMethod
    {
        Lines,
        Tubes,
        Ribbons
    };
    enum IntegrationDirection
    {
        Forward,
        Backward,
        Both
    };
    enum TerminationType
    {
        Distance,
        Time,
        Step
    };
    enum StreamlineAlgorithmType
    {
        LoadOnDemand,
        ParallelStaticDomains,
        MasterSlave
    };
    enum IntegrationType
    {
        DormandPrince,
        AdamsBashforth
    };
    enum OpacityType
    {
        None,
        Constant,
        Ramp,
        VariableRange
    };
    enum DisplayQuality
    {
        Low,
        Medium,
        High,
        Super
    };

    // These constructors are for objects of this class
    StreamlineAttributes();
    StreamlineAttributes(const StreamlineAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    StreamlineAttributes(private_tmfs_t tmfs);
    StreamlineAttributes(const StreamlineAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~StreamlineAttributes();

    virtual StreamlineAttributes& operator = (const StreamlineAttributes &obj);
    virtual bool operator == (const StreamlineAttributes &obj) const;
    virtual bool operator != (const StreamlineAttributes &obj) const;
private:
    void Init();
    void Copy(const StreamlineAttributes &obj);
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
    void SelectColorTableName();
    void SelectSingleColor();
    void SelectColoringVariable();
    void SelectOpacityVariable();

    // Property setting methods
    void SetSourceType(SourceType sourceType_);
    void SetMaxStepLength(double maxStepLength_);
    void SetTermination(double termination_);
    void SetPointSource(const double *pointSource_);
    void SetLineStart(const double *lineStart_);
    void SetLineEnd(const double *lineEnd_);
    void SetPlaneOrigin(const double *planeOrigin_);
    void SetPlaneNormal(const double *planeNormal_);
    void SetPlaneUpAxis(const double *planeUpAxis_);
    void SetPlaneRadius(double planeRadius_);
    void SetSphereOrigin(const double *sphereOrigin_);
    void SetSphereRadius(double sphereRadius_);
    void SetBoxExtents(const double *boxExtents_);
    void SetUseWholeBox(bool useWholeBox_);
    void SetPointList(const doubleVector &pointList_);
    void SetPointDensity(int pointDensity_);
    void SetDisplayMethod(DisplayMethod displayMethod_);
    void SetShowSeeds(bool showSeeds_);
    void SetShowHeads(bool showHeads_);
    void SetTubeRadius(double tubeRadius_);
    void SetRibbonWidth(double ribbonWidth_);
    void SetLineWidth(int lineWidth_);
    void SetColoringMethod(ColoringMethod coloringMethod_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetLegendFlag(bool legendFlag_);
    void SetLightingFlag(bool lightingFlag_);
    void SetStreamlineDirection(IntegrationDirection streamlineDirection_);
    void SetRelTol(double relTol_);
    void SetAbsTol(double absTol_);
    void SetTerminationType(TerminationType terminationType_);
    void SetIntegrationType(IntegrationType integrationType_);
    void SetStreamlineAlgorithmType(StreamlineAlgorithmType streamlineAlgorithmType_);
    void SetMaxStreamlineProcessCount(int maxStreamlineProcessCount_);
    void SetMaxDomainCacheSize(int maxDomainCacheSize_);
    void SetWorkGroupSize(int workGroupSize_);
    void SetPathlines(bool pathlines_);
    void SetColoringVariable(const std::string &coloringVariable_);
    void SetLegendMinFlag(bool legendMinFlag_);
    void SetLegendMaxFlag(bool legendMaxFlag_);
    void SetLegendMin(double legendMin_);
    void SetLegendMax(double legendMax_);
    void SetDisplayBegin(double displayBegin_);
    void SetDisplayEnd(double displayEnd_);
    void SetDisplayBeginFlag(bool displayBeginFlag_);
    void SetDisplayEndFlag(bool displayEndFlag_);
    void SetSeedDisplayRadius(double seedDisplayRadius_);
    void SetHeadDisplayRadius(double headDisplayRadius_);
    void SetOpacityType(OpacityType opacityType_);
    void SetOpacityVariable(const std::string &opacityVariable_);
    void SetOpacity(double opacity_);
    void SetOpacityVarMin(double opacityVarMin_);
    void SetOpacityVarMax(double opacityVarMax_);
    void SetOpacityVarMinFlag(bool opacityVarMinFlag_);
    void SetOpacityVarMaxFlag(bool opacityVarMaxFlag_);
    void SetTubeDisplayDensity(int tubeDisplayDensity_);
    void SetGeomDisplayQuality(DisplayQuality geomDisplayQuality_);

    // Property getting methods
    SourceType           GetSourceType() const;
    double               GetMaxStepLength() const;
    double               GetTermination() const;
    const double         *GetPointSource() const;
          double         *GetPointSource();
    const double         *GetLineStart() const;
          double         *GetLineStart();
    const double         *GetLineEnd() const;
          double         *GetLineEnd();
    const double         *GetPlaneOrigin() const;
          double         *GetPlaneOrigin();
    const double         *GetPlaneNormal() const;
          double         *GetPlaneNormal();
    const double         *GetPlaneUpAxis() const;
          double         *GetPlaneUpAxis();
    double               GetPlaneRadius() const;
    const double         *GetSphereOrigin() const;
          double         *GetSphereOrigin();
    double               GetSphereRadius() const;
    const double         *GetBoxExtents() const;
          double         *GetBoxExtents();
    bool                 GetUseWholeBox() const;
    const doubleVector   &GetPointList() const;
          doubleVector   &GetPointList();
    int                  GetPointDensity() const;
    DisplayMethod        GetDisplayMethod() const;
    bool                 GetShowSeeds() const;
    bool                 GetShowHeads() const;
    double               GetTubeRadius() const;
    double               GetRibbonWidth() const;
    int                  GetLineWidth() const;
    ColoringMethod       GetColoringMethod() const;
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    const ColorAttribute &GetSingleColor() const;
          ColorAttribute &GetSingleColor();
    bool                 GetLegendFlag() const;
    bool                 GetLightingFlag() const;
    IntegrationDirection GetStreamlineDirection() const;
    double               GetRelTol() const;
    double               GetAbsTol() const;
    TerminationType      GetTerminationType() const;
    IntegrationType      GetIntegrationType() const;
    StreamlineAlgorithmType GetStreamlineAlgorithmType() const;
    int                  GetMaxStreamlineProcessCount() const;
    int                  GetMaxDomainCacheSize() const;
    int                  GetWorkGroupSize() const;
    bool                 GetPathlines() const;
    const std::string    &GetColoringVariable() const;
          std::string    &GetColoringVariable();
    bool                 GetLegendMinFlag() const;
    bool                 GetLegendMaxFlag() const;
    double               GetLegendMin() const;
    double               GetLegendMax() const;
    double               GetDisplayBegin() const;
    double               GetDisplayEnd() const;
    bool                 GetDisplayBeginFlag() const;
    bool                 GetDisplayEndFlag() const;
    double               GetSeedDisplayRadius() const;
    double               GetHeadDisplayRadius() const;
    OpacityType          GetOpacityType() const;
    const std::string    &GetOpacityVariable() const;
          std::string    &GetOpacityVariable();
    double               GetOpacity() const;
    double               GetOpacityVarMin() const;
    double               GetOpacityVarMax() const;
    bool                 GetOpacityVarMinFlag() const;
    bool                 GetOpacityVarMaxFlag() const;
    int                  GetTubeDisplayDensity() const;
    DisplayQuality       GetGeomDisplayQuality() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string SourceType_ToString(SourceType);
    static bool SourceType_FromString(const std::string &, SourceType &);
protected:
    static std::string SourceType_ToString(int);
public:
    static std::string ColoringMethod_ToString(ColoringMethod);
    static bool ColoringMethod_FromString(const std::string &, ColoringMethod &);
protected:
    static std::string ColoringMethod_ToString(int);
public:
    static std::string DisplayMethod_ToString(DisplayMethod);
    static bool DisplayMethod_FromString(const std::string &, DisplayMethod &);
protected:
    static std::string DisplayMethod_ToString(int);
public:
    static std::string IntegrationDirection_ToString(IntegrationDirection);
    static bool IntegrationDirection_FromString(const std::string &, IntegrationDirection &);
protected:
    static std::string IntegrationDirection_ToString(int);
public:
    static std::string TerminationType_ToString(TerminationType);
    static bool TerminationType_FromString(const std::string &, TerminationType &);
protected:
    static std::string TerminationType_ToString(int);
public:
    static std::string StreamlineAlgorithmType_ToString(StreamlineAlgorithmType);
    static bool StreamlineAlgorithmType_FromString(const std::string &, StreamlineAlgorithmType &);
protected:
    static std::string StreamlineAlgorithmType_ToString(int);
public:
    static std::string IntegrationType_ToString(IntegrationType);
    static bool IntegrationType_FromString(const std::string &, IntegrationType &);
protected:
    static std::string IntegrationType_ToString(int);
public:
    static std::string OpacityType_ToString(OpacityType);
    static bool OpacityType_FromString(const std::string &, OpacityType &);
protected:
    static std::string OpacityType_ToString(int);
public:
    static std::string DisplayQuality_ToString(DisplayQuality);
    static bool DisplayQuality_FromString(const std::string &, DisplayQuality &);
protected:
    static std::string DisplayQuality_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const StreamlineAttributes &) const;
    virtual void  ProcessOldVersions(DataNode *parentNode, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_sourceType = 0,
        ID_maxStepLength,
        ID_termination,
        ID_pointSource,
        ID_lineStart,
        ID_lineEnd,
        ID_planeOrigin,
        ID_planeNormal,
        ID_planeUpAxis,
        ID_planeRadius,
        ID_sphereOrigin,
        ID_sphereRadius,
        ID_boxExtents,
        ID_useWholeBox,
        ID_pointList,
        ID_pointDensity,
        ID_displayMethod,
        ID_showSeeds,
        ID_showHeads,
        ID_tubeRadius,
        ID_ribbonWidth,
        ID_lineWidth,
        ID_coloringMethod,
        ID_colorTableName,
        ID_singleColor,
        ID_legendFlag,
        ID_lightingFlag,
        ID_streamlineDirection,
        ID_relTol,
        ID_absTol,
        ID_terminationType,
        ID_integrationType,
        ID_streamlineAlgorithmType,
        ID_maxStreamlineProcessCount,
        ID_maxDomainCacheSize,
        ID_workGroupSize,
        ID_pathlines,
        ID_coloringVariable,
        ID_legendMinFlag,
        ID_legendMaxFlag,
        ID_legendMin,
        ID_legendMax,
        ID_displayBegin,
        ID_displayEnd,
        ID_displayBeginFlag,
        ID_displayEndFlag,
        ID_seedDisplayRadius,
        ID_headDisplayRadius,
        ID_opacityType,
        ID_opacityVariable,
        ID_opacity,
        ID_opacityVarMin,
        ID_opacityVarMax,
        ID_opacityVarMinFlag,
        ID_opacityVarMaxFlag,
        ID_tubeDisplayDensity,
        ID_geomDisplayQuality,
        ID__LAST
    };

private:
    int            sourceType;
    double         maxStepLength;
    double         termination;
    double         pointSource[3];
    double         lineStart[3];
    double         lineEnd[3];
    double         planeOrigin[3];
    double         planeNormal[3];
    double         planeUpAxis[3];
    double         planeRadius;
    double         sphereOrigin[3];
    double         sphereRadius;
    double         boxExtents[6];
    bool           useWholeBox;
    doubleVector   pointList;
    int            pointDensity;
    int            displayMethod;
    bool           showSeeds;
    bool           showHeads;
    double         tubeRadius;
    double         ribbonWidth;
    int            lineWidth;
    int            coloringMethod;
    std::string    colorTableName;
    ColorAttribute singleColor;
    bool           legendFlag;
    bool           lightingFlag;
    int            streamlineDirection;
    double         relTol;
    double         absTol;
    int            terminationType;
    int            integrationType;
    int            streamlineAlgorithmType;
    int            maxStreamlineProcessCount;
    int            maxDomainCacheSize;
    int            workGroupSize;
    bool           pathlines;
    std::string    coloringVariable;
    bool           legendMinFlag;
    bool           legendMaxFlag;
    double         legendMin;
    double         legendMax;
    double         displayBegin;
    double         displayEnd;
    bool           displayBeginFlag;
    bool           displayEndFlag;
    double         seedDisplayRadius;
    double         headDisplayRadius;
    int            opacityType;
    std::string    opacityVariable;
    double         opacity;
    double         opacityVarMin;
    double         opacityVarMax;
    bool           opacityVarMinFlag;
    bool           opacityVarMaxFlag;
    int            tubeDisplayDensity;
    int            geomDisplayQuality;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define STREAMLINEATTRIBUTES_TMFS "iddDDDDDDdDdDbd*iibbddiisabbiddiiiiiibsbbddddbbddisdddbbii"

#endif
