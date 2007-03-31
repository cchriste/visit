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
// Creation:   Thu Oct 9 13:31:00 PST 2003
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

    StreamlineAttributes();
    StreamlineAttributes(const StreamlineAttributes &obj);
    virtual ~StreamlineAttributes();

    virtual void operator = (const StreamlineAttributes &obj);
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
    void SetStepLength(double stepLength_);
    void SetMaxTime(double maxTime_);
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
    void SetPointDensity(int pointDensity_);
    void SetShowTube(bool showTube_);
    void SetShowStart(bool showStart_);
    void SetTubeRadius(double tubeRadius_);
    void SetLineWidth(int lineWidth_);
    void SetColorBySpeed(bool colorBySpeed_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetSingleColor(const ColorAttribute &singleColor_);
    void SetLegendFlag(bool legendFlag_);
    void SetLightingFlag(bool lightingFlag_);

    // Property getting methods
    SourceType           GetSourceType() const;
    double               GetStepLength() const;
    double               GetMaxTime() const;
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
    int                  GetPointDensity() const;
    bool                 GetShowTube() const;
    bool                 GetShowStart() const;
    double               GetTubeRadius() const;
    int                  GetLineWidth() const;
    bool                 GetColorBySpeed() const;
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    const ColorAttribute &GetSingleColor() const;
          ColorAttribute &GetSingleColor();
    bool                 GetLegendFlag() const;
    bool                 GetLightingFlag() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string SourceType_ToString(SourceType);
    static bool SourceType_FromString(const std::string &, SourceType &);
protected:
    static std::string SourceType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const StreamlineAttributes &) const;
private:
    int            sourceType;
    double         stepLength;
    double         maxTime;
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
    int            pointDensity;
    bool           showTube;
    bool           showStart;
    double         tubeRadius;
    int            lineWidth;
    bool           colorBySpeed;
    std::string    colorTableName;
    ColorAttribute singleColor;
    bool           legendFlag;
    bool           lightingFlag;
};

#endif
