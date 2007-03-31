#ifndef MESHATTRIBUTES_H
#define MESHATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: MeshAttributes
//
// Purpose:
//    Attributes for the mesh plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Aug 21 23:06:49 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class MeshAttributes : public AttributeSubject
{
public:
    enum PointType
    {
        Box,
        Axis,
        Icosahedron
    };

    MeshAttributes();
    MeshAttributes(const MeshAttributes &obj);
    virtual ~MeshAttributes();

    virtual void operator = (const MeshAttributes &obj);
    virtual bool operator == (const MeshAttributes &obj) const;
    virtual bool operator != (const MeshAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectMeshColor();
    void SelectOpaqueColor();
    void SelectPointSizeVar();

    // Property setting methods
    void SetLegendFlag(bool legendFlag_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetMeshColor(const ColorAttribute &meshColor_);
    void SetOutlineOnlyFlag(bool outlineOnlyFlag_);
    void SetErrorTolerance(double errorTolerance_);
    void SetOpaqueFlag(bool opaqueFlag_);
    void SetPointSize(double pointSize_);
    void SetOpaqueColor(const ColorAttribute &opaqueColor_);
    void SetBackgroundFlag(bool backgroundFlag_);
    void SetForegroundFlag(bool foregroundFlag_);
    void SetSmoothingLevel(int smoothingLevel_);
    void SetPointSizeVarEnabled(bool pointSizeVarEnabled_);
    void SetPointSizeVar(const std::string &pointSizeVar_);
    void SetPointType(PointType pointType_);

    // Property getting methods
    bool                 GetLegendFlag() const;
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    const ColorAttribute &GetMeshColor() const;
          ColorAttribute &GetMeshColor();
    bool                 GetOutlineOnlyFlag() const;
    double               GetErrorTolerance() const;
    bool                 GetOpaqueFlag() const;
    double               GetPointSize() const;
    const ColorAttribute &GetOpaqueColor() const;
          ColorAttribute &GetOpaqueColor();
    bool                 GetBackgroundFlag() const;
    bool                 GetForegroundFlag() const;
    int                  GetSmoothingLevel() const;
    bool                 GetPointSizeVarEnabled() const;
    const std::string    &GetPointSizeVar() const;
          std::string    &GetPointSizeVar();
    PointType            GetPointType() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string PointType_ToString(PointType);
    static bool PointType_FromString(const std::string &, PointType &);
protected:
    static std::string PointType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const MeshAttributes &obj);
private:
    bool           legendFlag;
    int            lineStyle;
    int            lineWidth;
    ColorAttribute meshColor;
    bool           outlineOnlyFlag;
    double         errorTolerance;
    bool           opaqueFlag;
    double         pointSize;
    ColorAttribute opaqueColor;
    bool           backgroundFlag;
    bool           foregroundFlag;
    int            smoothingLevel;
    bool           pointSizeVarEnabled;
    std::string    pointSizeVar;
    int            pointType;
};

#endif
