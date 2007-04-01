#ifndef VECTORATTRIBUTES_H
#define VECTORATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: VectorAttributes
//
// Purpose:
//    Attributes for the vector plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Jun 21 15:36:18 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class VectorAttributes : public AttributeSubject
{
public:
    enum OriginType
    {
        Head,
        Middle,
        Tail
    };

    VectorAttributes();
    VectorAttributes(const VectorAttributes &obj);
    virtual ~VectorAttributes();

    virtual void operator = (const VectorAttributes &obj);
    virtual bool operator == (const VectorAttributes &obj) const;
    virtual bool operator != (const VectorAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVectorColor();
    void SelectColorTableName();

    // Property setting methods
    void SetUseStride(bool useStride_);
    void SetStride(int stride_);
    void SetNVectors(int nVectors_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetScale(double scale_);
    void SetHeadSize(double headSize_);
    void SetHeadOn(bool headOn_);
    void SetColorByMag(bool colorByMag_);
    void SetUseLegend(bool useLegend_);
    void SetVectorColor(const ColorAttribute &vectorColor_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetVectorOrigin(OriginType vectorOrigin_);

    // Property getting methods
    bool                 GetUseStride() const;
    int                  GetStride() const;
    int                  GetNVectors() const;
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    double               GetScale() const;
    double               GetHeadSize() const;
    bool                 GetHeadOn() const;
    bool                 GetColorByMag() const;
    bool                 GetUseLegend() const;
    const ColorAttribute &GetVectorColor() const;
          ColorAttribute &GetVectorColor();
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    OriginType           GetVectorOrigin() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string OriginType_ToString(OriginType);
    static bool OriginType_FromString(const std::string &, OriginType &);
protected:
    static std::string OriginType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const VectorAttributes &obj);
private:
    bool           useStride;
    int            stride;
    int            nVectors;
    int            lineStyle;
    int            lineWidth;
    double         scale;
    double         headSize;
    bool           headOn;
    bool           colorByMag;
    bool           useLegend;
    ColorAttribute vectorColor;
    std::string    colorTableName;
    int            vectorOrigin;
};

#endif
