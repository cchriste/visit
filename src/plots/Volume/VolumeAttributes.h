#ifndef VOLUMEATTRIBUTES_H
#define VOLUMEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorControlPointList.h>
#include <GaussianControlPointList.h>

// ****************************************************************************
// Class: VolumeAttributes
//
// Purpose:
//    This class contains the plot attributes for the volume plot.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Oct 2 12:43:07 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class VolumeAttributes : public AttributeSubject
{
public:
    enum Renderer
    {
        Splatting,
        Texture3D
    };
    enum GradientType
    {
        CenteredDifferences,
        SobelOperator
    };

    VolumeAttributes();
    VolumeAttributes(const VolumeAttributes &obj);
    virtual ~VolumeAttributes();

    virtual void operator = (const VolumeAttributes &obj);
    virtual bool operator == (const VolumeAttributes &obj) const;
    virtual bool operator != (const VolumeAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColorControlPoints();
    void SelectOpacityControlPoints();
    void SelectOpacityVariable();
    void SelectFreeformOpacity();

    // Property setting methods
    void SetLegendFlag(bool legendFlag_);
    void SetLightingFlag(bool lightingFlag_);
    void SetColorControlPoints(const ColorControlPointList &colorControlPoints_);
    void SetOpacityAttenuation(float opacityAttenuation_);
    void SetFreeformFlag(bool freeformFlag_);
    void SetOpacityControlPoints(const GaussianControlPointList &opacityControlPoints_);
    void SetResampleTarget(int resampleTarget_);
    void SetOpacityVariable(const std::string &opacityVariable_);
    void SetFreeformOpacity(const unsigned char *freeformOpacity_);
    void SetDoSoftware(bool doSoftware_);
    void SetUseColorVarMin(bool useColorVarMin_);
    void SetColorVarMin(float colorVarMin_);
    void SetUseColorVarMax(bool useColorVarMax_);
    void SetColorVarMax(float colorVarMax_);
    void SetUseOpacityVarMin(bool useOpacityVarMin_);
    void SetOpacityVarMin(float opacityVarMin_);
    void SetUseOpacityVarMax(bool useOpacityVarMax_);
    void SetOpacityVarMax(float opacityVarMax_);
    void SetSmoothData(bool smoothData_);
    void SetSamplesPerRay(int samplesPerRay_);
    void SetRendererType(Renderer rendererType_);
    void SetGradientType(GradientType gradientType_);
    void SetNum3DSlices(int num3DSlices_);

    // Property getting methods
    bool                           GetLegendFlag() const;
    bool                           GetLightingFlag() const;
    const ColorControlPointList    &GetColorControlPoints() const;
          ColorControlPointList    &GetColorControlPoints();
    float                          GetOpacityAttenuation() const;
    bool                           GetFreeformFlag() const;
    const GaussianControlPointList &GetOpacityControlPoints() const;
          GaussianControlPointList &GetOpacityControlPoints();
    int                            GetResampleTarget() const;
    const std::string              &GetOpacityVariable() const;
          std::string              &GetOpacityVariable();
    const unsigned char            *GetFreeformOpacity() const;
          unsigned char            *GetFreeformOpacity();
    bool                           GetDoSoftware() const;
    bool                           GetUseColorVarMin() const;
    float                          GetColorVarMin() const;
    bool                           GetUseColorVarMax() const;
    float                          GetColorVarMax() const;
    bool                           GetUseOpacityVarMin() const;
    float                          GetOpacityVarMin() const;
    bool                           GetUseOpacityVarMax() const;
    float                          GetOpacityVarMax() const;
    bool                           GetSmoothData() const;
    int                            GetSamplesPerRay() const;
    Renderer                       GetRendererType() const;
    GradientType                   GetGradientType() const;
    int                            GetNum3DSlices() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Renderer_ToString(Renderer);
    static bool Renderer_FromString(const std::string &, Renderer &);
protected:
    static std::string Renderer_ToString(int);
public:
    static std::string GradientType_ToString(GradientType);
    static bool GradientType_FromString(const std::string &, GradientType &);
protected:
    static std::string GradientType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const VolumeAttributes &obj) const;
    bool GradientWontChange(const VolumeAttributes &obj) const;
    void GetTransferFunction(unsigned char *rgba) const;
    void SetDefaultColorControlPoints();
    void GetGaussianOpacities(unsigned char *alphas) const;
    void GetOpacities(unsigned char *alphas);
    void SetSmoothingFlag(bool val);
    bool GetSmoothingFlag() const;
    void SetEqualSpacingFlag(bool val);
    bool GetEqualSpacingFlag() const;
private:
    bool                     legendFlag;
    bool                     lightingFlag;
    ColorControlPointList    colorControlPoints;
    float                    opacityAttenuation;
    bool                     freeformFlag;
    GaussianControlPointList opacityControlPoints;
    int                      resampleTarget;
    std::string              opacityVariable;
    unsigned char            freeformOpacity[256];
    bool                     doSoftware;
    bool                     useColorVarMin;
    float                    colorVarMin;
    bool                     useColorVarMax;
    float                    colorVarMax;
    bool                     useOpacityVarMin;
    float                    opacityVarMin;
    bool                     useOpacityVarMax;
    float                    opacityVarMax;
    bool                     smoothData;
    int                      samplesPerRay;
    int                      rendererType;
    int                      gradientType;
    int                      num3DSlices;
};

#endif
