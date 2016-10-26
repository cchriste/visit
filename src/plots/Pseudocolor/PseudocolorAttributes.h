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

#ifndef PSEUDOCOLORATTRIBUTES_H
#define PSEUDOCOLORATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>
#include <visitstream.h>

// ****************************************************************************
// Class: PseudocolorAttributes
//
// Purpose:
//    Attributes for the pseudocolor plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class PseudocolorAttributes : public AttributeSubject
{
public:
    enum Scaling
    {
        Linear,
        Log,
        Skew
    };
    enum LimitsMode
    {
        OriginalData,
        CurrentPlot
    };
    enum Centering
    {
        Natural,
        Nodal,
        Zonal
    };
    enum OpacityType
    {
        ColorTable,
        FullyOpaque,
        Constant,
        Ramp,
        VariableRange
    };
    enum PointType
    {
        Box,
        Axis,
        Icosahedron,
        Octahedron,
        Tetrahedron,
        SphereGeometry,
        Point,
        Sphere
    };
    enum LineType
    {
        Line,
        Tube,
        Ribbon
    };
    enum EndPointStyle
    {
        EndPointNone,
        EndPointSphere,
        EndPointCone
    };
    enum SizeType
    {
        Absolute,
        FractionOfBBox
    };

    // These constructors are for objects of this class
    PseudocolorAttributes();
    PseudocolorAttributes(const PseudocolorAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    PseudocolorAttributes(private_tmfs_t tmfs);
    PseudocolorAttributes(const PseudocolorAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~PseudocolorAttributes();

    virtual PseudocolorAttributes& operator = (const PseudocolorAttributes &obj);
    virtual bool operator == (const PseudocolorAttributes &obj) const;
    virtual bool operator != (const PseudocolorAttributes &obj) const;
private:
    void Init();
    void Copy(const PseudocolorAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColorTableName();
    void SelectOpacityVariable();
    void SelectPointSizeVar();
    void SelectTubeRadiusVar();
    void SelectEndPointRadiusVar();
    void SelectWireframeColor();
    void SelectPointColor();

    // Property setting methods
    void SetScaling(Scaling scaling_);
    void SetSkewFactor(double skewFactor_);
    void SetLimitsMode(LimitsMode limitsMode_);
    void SetMinFlag(bool minFlag_);
    void SetMin(double min_);
    void SetMaxFlag(bool maxFlag_);
    void SetMax(double max_);
    void SetCentering(Centering centering_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetInvertColorTable(bool invertColorTable_);
    void SetOpacityType(OpacityType opacityType_);
    void SetOpacityVariable(const std::string &opacityVariable_);
    void SetOpacity(double opacity_);
    void SetOpacityVarMin(double opacityVarMin_);
    void SetOpacityVarMax(double opacityVarMax_);
    void SetOpacityVarMinFlag(bool opacityVarMinFlag_);
    void SetOpacityVarMaxFlag(bool opacityVarMaxFlag_);
    void SetPointSize(double pointSize_);
    void SetPointType(PointType pointType_);
    void SetPointSizeVarEnabled(bool pointSizeVarEnabled_);
    void SetPointSizeVar(const std::string &pointSizeVar_);
    void SetPointSizePixels(int pointSizePixels_);
    void SetLineStyle(int lineStyle_);
    void SetLineType(LineType lineType_);
    void SetLineWidth(int lineWidth_);
    void SetTubeResolution(int tubeResolution_);
    void SetTubeRadiusSizeType(SizeType tubeRadiusSizeType_);
    void SetTubeRadiusAbsolute(double tubeRadiusAbsolute_);
    void SetTubeRadiusBBox(double tubeRadiusBBox_);
    void SetTubeRadiusVarEnabled(bool tubeRadiusVarEnabled_);
    void SetTubeRadiusVar(const std::string &tubeRadiusVar_);
    void SetTubeRadiusVarRatio(double tubeRadiusVarRatio_);
    void SetTailStyle(EndPointStyle tailStyle_);
    void SetHeadStyle(EndPointStyle headStyle_);
    void SetEndPointRadiusSizeType(SizeType endPointRadiusSizeType_);
    void SetEndPointRadiusAbsolute(double endPointRadiusAbsolute_);
    void SetEndPointRadiusBBox(double endPointRadiusBBox_);
    void SetEndPointResolution(int endPointResolution_);
    void SetEndPointRatio(double endPointRatio_);
    void SetEndPointRadiusVarEnabled(bool endPointRadiusVarEnabled_);
    void SetEndPointRadiusVar(const std::string &endPointRadiusVar_);
    void SetEndPointRadiusVarRatio(double endPointRadiusVarRatio_);
    void SetRenderSurfaces(int renderSurfaces_);
    void SetRenderWireframe(int renderWireframe_);
    void SetRenderPoints(int renderPoints_);
    void SetSmoothingLevel(int smoothingLevel_);
    void SetLegendFlag(bool legendFlag_);
    void SetLightingFlag(bool lightingFlag_);
    void SetWireframeColor(const ColorAttribute &wireframeColor_);
    void SetPointColor(const ColorAttribute &pointColor_);

    // Property getting methods
    Scaling              GetScaling() const;
    double               GetSkewFactor() const;
    LimitsMode           GetLimitsMode() const;
    bool                 GetMinFlag() const;
    double               GetMin() const;
    bool                 GetMaxFlag() const;
    double               GetMax() const;
    Centering            GetCentering() const;
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    bool                 GetInvertColorTable() const;
    OpacityType          GetOpacityType() const;
    const std::string    &GetOpacityVariable() const;
          std::string    &GetOpacityVariable();
    double               GetOpacity() const;
    double               GetOpacityVarMin() const;
    double               GetOpacityVarMax() const;
    bool                 GetOpacityVarMinFlag() const;
    bool                 GetOpacityVarMaxFlag() const;
    double               GetPointSize() const;
    PointType            GetPointType() const;
    bool                 GetPointSizeVarEnabled() const;
    const std::string    &GetPointSizeVar() const;
          std::string    &GetPointSizeVar();
    int                  GetPointSizePixels() const;
    int                  GetLineStyle() const;
    LineType             GetLineType() const;
    int                  GetLineWidth() const;
    int                  GetTubeResolution() const;
    SizeType             GetTubeRadiusSizeType() const;
    double               GetTubeRadiusAbsolute() const;
    double               GetTubeRadiusBBox() const;
    bool                 GetTubeRadiusVarEnabled() const;
    const std::string    &GetTubeRadiusVar() const;
          std::string    &GetTubeRadiusVar();
    double               GetTubeRadiusVarRatio() const;
    EndPointStyle        GetTailStyle() const;
    EndPointStyle        GetHeadStyle() const;
    SizeType             GetEndPointRadiusSizeType() const;
    double               GetEndPointRadiusAbsolute() const;
    double               GetEndPointRadiusBBox() const;
    int                  GetEndPointResolution() const;
    double               GetEndPointRatio() const;
    bool                 GetEndPointRadiusVarEnabled() const;
    const std::string    &GetEndPointRadiusVar() const;
          std::string    &GetEndPointRadiusVar();
    double               GetEndPointRadiusVarRatio() const;
    int                  GetRenderSurfaces() const;
    int                  GetRenderWireframe() const;
    int                  GetRenderPoints() const;
    int                  GetSmoothingLevel() const;
    bool                 GetLegendFlag() const;
    bool                 GetLightingFlag() const;
    const ColorAttribute &GetWireframeColor() const;
          ColorAttribute &GetWireframeColor();
    const ColorAttribute &GetPointColor() const;
          ColorAttribute &GetPointColor();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Scaling_ToString(Scaling);
    static bool Scaling_FromString(const std::string &, Scaling &);
protected:
    static std::string Scaling_ToString(int);
public:
    static std::string LimitsMode_ToString(LimitsMode);
    static bool LimitsMode_FromString(const std::string &, LimitsMode &);
protected:
    static std::string LimitsMode_ToString(int);
public:
    static std::string Centering_ToString(Centering);
    static bool Centering_FromString(const std::string &, Centering &);
protected:
    static std::string Centering_ToString(int);
public:
    static std::string OpacityType_ToString(OpacityType);
    static bool OpacityType_FromString(const std::string &, OpacityType &);
protected:
    static std::string OpacityType_ToString(int);
public:
    static std::string PointType_ToString(PointType);
    static bool PointType_FromString(const std::string &, PointType &);
protected:
    static std::string PointType_ToString(int);
public:
    static std::string LineType_ToString(LineType);
    static bool LineType_FromString(const std::string &, LineType &);
protected:
    static std::string LineType_ToString(int);
public:
    static std::string EndPointStyle_ToString(EndPointStyle);
    static bool EndPointStyle_FromString(const std::string &, EndPointStyle &);
protected:
    static std::string EndPointStyle_ToString(int);
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
    bool ChangesRequireRecalculation(const PseudocolorAttributes &) const;
    void Print(ostream &, bool) const;
    virtual void ProcessOldVersions(DataNode *parentNode, const char *configVersion);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_scaling = 0,
        ID_skewFactor,
        ID_limitsMode,
        ID_minFlag,
        ID_min,
        ID_maxFlag,
        ID_max,
        ID_centering,
        ID_colorTableName,
        ID_invertColorTable,
        ID_opacityType,
        ID_opacityVariable,
        ID_opacity,
        ID_opacityVarMin,
        ID_opacityVarMax,
        ID_opacityVarMinFlag,
        ID_opacityVarMaxFlag,
        ID_pointSize,
        ID_pointType,
        ID_pointSizeVarEnabled,
        ID_pointSizeVar,
        ID_pointSizePixels,
        ID_lineStyle,
        ID_lineType,
        ID_lineWidth,
        ID_tubeResolution,
        ID_tubeRadiusSizeType,
        ID_tubeRadiusAbsolute,
        ID_tubeRadiusBBox,
        ID_tubeRadiusVarEnabled,
        ID_tubeRadiusVar,
        ID_tubeRadiusVarRatio,
        ID_tailStyle,
        ID_headStyle,
        ID_endPointRadiusSizeType,
        ID_endPointRadiusAbsolute,
        ID_endPointRadiusBBox,
        ID_endPointResolution,
        ID_endPointRatio,
        ID_endPointRadiusVarEnabled,
        ID_endPointRadiusVar,
        ID_endPointRadiusVarRatio,
        ID_renderSurfaces,
        ID_renderWireframe,
        ID_renderPoints,
        ID_smoothingLevel,
        ID_legendFlag,
        ID_lightingFlag,
        ID_wireframeColor,
        ID_pointColor,
        ID__LAST
    };

private:
    int            scaling;
    double         skewFactor;
    int            limitsMode;
    bool           minFlag;
    double         min;
    bool           maxFlag;
    double         max;
    int            centering;
    std::string    colorTableName;
    bool           invertColorTable;
    int            opacityType;
    std::string    opacityVariable;
    double         opacity;
    double         opacityVarMin;
    double         opacityVarMax;
    bool           opacityVarMinFlag;
    bool           opacityVarMaxFlag;
    double         pointSize;
    int            pointType;
    bool           pointSizeVarEnabled;
    std::string    pointSizeVar;
    int            pointSizePixels;
    int            lineStyle;
    int            lineType;
    int            lineWidth;
    int            tubeResolution;
    int            tubeRadiusSizeType;
    double         tubeRadiusAbsolute;
    double         tubeRadiusBBox;
    bool           tubeRadiusVarEnabled;
    std::string    tubeRadiusVar;
    double         tubeRadiusVarRatio;
    int            tailStyle;
    int            headStyle;
    int            endPointRadiusSizeType;
    double         endPointRadiusAbsolute;
    double         endPointRadiusBBox;
    int            endPointResolution;
    double         endPointRatio;
    bool           endPointRadiusVarEnabled;
    std::string    endPointRadiusVar;
    double         endPointRadiusVarRatio;
    int            renderSurfaces;
    int            renderWireframe;
    int            renderPoints;
    int            smoothingLevel;
    bool           legendFlag;
    bool           lightingFlag;
    ColorAttribute wireframeColor;
    ColorAttribute pointColor;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define PSEUDOCOLORATTRIBUTES_TMFS "idibdbdisbisdddbbdibsiiiiiiddbsdiiiddidbsdiiiibbaa"

#endif
