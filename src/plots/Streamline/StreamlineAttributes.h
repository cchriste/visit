/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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
        SpecifiedBox
    };
    enum ColoringMethod
    {
        Solid,
        ColorBySpeed,
        ColorByVorticity,
        ColorByLength,
        ColorByTime,
        ColorBySeedPointID
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
    enum IntegrationType
    {
        DormandPrince,
        AdamsBashforth
    };
    enum StreamlineAlgorithmType
    {
        LoadOnDemand,
        ParallelStaticDomains,
        MasterSlave
    };

    StreamlineAttributes();
    StreamlineAttributes(const StreamlineAttributes &obj);
    virtual ~StreamlineAttributes();

    virtual StreamlineAttributes& operator = (const StreamlineAttributes &obj);
    virtual bool operator == (const StreamlineAttributes &obj) const;
    virtual bool operator != (const StreamlineAttributes &obj) const;

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
    void SelectColorTableName();
    void SelectSingleColor();

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
    void SetPointDensity(int pointDensity_);
    void SetDisplayMethod(DisplayMethod displayMethod_);
    void SetShowStart(bool showStart_);
    void SetRadius(double radius_);
    void SetLineWidth(int lineWidth_);
    void SetColoringMethod(ColoringMethod coloringMethod_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetLegendFlag(bool legendFlag_);
    void SetLightingFlag(bool lightingFlag_);
    void SetStreamlineDirection(IntegrationDirection StreamlineDirection_);
    void SetRelTol(double relTol_);
    void SetAbsTol(double absTol_);
    void SetTerminationType(TerminationType terminationType_);
    void SetIntegrationType(IntegrationType integrationType_);
    void SetStreamlineAlgorithmType(StreamlineAlgorithmType streamlineAlgorithmType_);
    void SetMaxStreamlineProcessCount(int maxStreamlineProcessCount_);
    void SetMaxDomainCacheSize(int maxDomainCacheSize_);
    void SetWorkGroupSize(int workGroupSize_);

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
    int                  GetPointDensity() const;
    DisplayMethod        GetDisplayMethod() const;
    bool                 GetShowStart() const;
    double               GetRadius() const;
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
    static std::string IntegrationType_ToString(IntegrationType);
    static bool IntegrationType_FromString(const std::string &, IntegrationType &);
protected:
    static std::string IntegrationType_ToString(int);
public:
    static std::string StreamlineAlgorithmType_ToString(StreamlineAlgorithmType);
    static bool StreamlineAlgorithmType_FromString(const std::string &, StreamlineAlgorithmType &);
protected:
    static std::string StreamlineAlgorithmType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const StreamlineAttributes &) const;

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
        ID_pointDensity,
        ID_displayMethod,
        ID_showStart,
        ID_radius,
        ID_lineWidth,
        ID_coloringMethod,
        ID_colorTableName,
        ID_singleColor,
        ID_legendFlag,
        ID_lightingFlag,
        ID_StreamlineDirection,
        ID_relTol,
        ID_absTol,
        ID_terminationType,
        ID_integrationType,
        ID_streamlineAlgorithmType,
        ID_maxStreamlineProcessCount,
        ID_maxDomainCacheSize,
        ID_workGroupSize
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
    int            pointDensity;
    int            displayMethod;
    bool           showStart;
    double         radius;
    int            lineWidth;
    int            coloringMethod;
    std::string    colorTableName;
    ColorAttribute singleColor;
    bool           legendFlag;
    bool           lightingFlag;
    int            StreamlineDirection;
    double         relTol;
    double         absTol;
    int            terminationType;
    int            integrationType;
    int            streamlineAlgorithmType;
    int            maxStreamlineProcessCount;
    int            maxDomainCacheSize;
    int            workGroupSize;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
};

#endif
