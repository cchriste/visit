#ifndef CURVEATTRIBUTES_H
#define CURVEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <ColorAttribute.h>
#include <visitstream.h>

// ****************************************************************************
// Class: CurveAttributes
//
// Purpose:
//    Attributes for the xy plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Dec 23 13:34:17 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class CurveAttributes : public AttributeSubject
{
public:
    CurveAttributes();
    CurveAttributes(const CurveAttributes &obj);
    virtual ~CurveAttributes();

    virtual void operator = (const CurveAttributes &obj);
    virtual bool operator == (const CurveAttributes &obj) const;
    virtual bool operator != (const CurveAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectColor();
    void SelectDesignator();

    // Property setting methods
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetColor(const ColorAttribute &color_);
    void SetShowLabels(bool showLabels_);
    void SetDesignator(const std::string &designator_);
    void SetShowPoints(bool showPoints_);
    void SetPointSize(double pointSize_);

    // Property getting methods
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    const ColorAttribute &GetColor() const;
          ColorAttribute &GetColor();
    bool                 GetShowLabels() const;
    const std::string    &GetDesignator() const;
          std::string    &GetDesignator();
    bool                 GetShowPoints() const;
    double               GetPointSize() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const CurveAttributes &) const;
    void Print(ostream &, bool) const;
private:
    int            lineStyle;
    int            lineWidth;
    ColorAttribute color;
    bool           showLabels;
    std::string    designator;
    bool           showPoints;
    double         pointSize;
};

#endif
