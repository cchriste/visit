#ifndef TENSORATTRIBUTES_H
#define TENSORATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>

// ****************************************************************************
// Class: TensorAttributes
//
// Purpose:
//    Attributes for the tensor plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Oct 9 13:31:08 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class TensorAttributes : public AttributeSubject
{
public:
    TensorAttributes();
    TensorAttributes(const TensorAttributes &obj);
    virtual ~TensorAttributes();

    virtual void operator = (const TensorAttributes &obj);
    virtual bool operator == (const TensorAttributes &obj) const;
    virtual bool operator != (const TensorAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTensorColor();
    void SelectColorTableName();

    // Property setting methods
    void SetUseStride(bool useStride_);
    void SetStride(int stride_);
    void SetNTensors(int nTensors_);
    void SetScale(double scale_);
    void SetColorByEigenvalues(bool colorByEigenvalues_);
    void SetUseLegend(bool useLegend_);
    void SetTensorColor(const ColorAttribute &tensorColor_);
    void SetColorTableName(const std::string &colorTableName_);

    // Property getting methods
    bool                 GetUseStride() const;
    int                  GetStride() const;
    int                  GetNTensors() const;
    double               GetScale() const;
    bool                 GetColorByEigenvalues() const;
    bool                 GetUseLegend() const;
    const ColorAttribute &GetTensorColor() const;
          ColorAttribute &GetTensorColor();
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const TensorAttributes &obj);
private:
    bool           useStride;
    int            stride;
    int            nTensors;
    double         scale;
    bool           colorByEigenvalues;
    bool           useLegend;
    ColorAttribute tensorColor;
    std::string    colorTableName;
};

#endif
