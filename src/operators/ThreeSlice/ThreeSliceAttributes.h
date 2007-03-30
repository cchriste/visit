#ifndef THREESLICEATTRIBUTES_H
#define THREESLICEATTRIBUTES_H
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ThreeSliceAttributes
//
// Purpose:
//    This class contains attributes for the threeslice operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Jun 9 13:18:30 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class ThreeSliceAttributes : public AttributeSubject
{
public:
    ThreeSliceAttributes();
    ThreeSliceAttributes(const ThreeSliceAttributes &obj);
    virtual ~ThreeSliceAttributes();

    virtual void operator = (const ThreeSliceAttributes &obj);
    virtual bool operator == (const ThreeSliceAttributes &obj);
    virtual bool operator != (const ThreeSliceAttributes &obj);

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetX(float x_);
    void SetY(float y_);
    void SetZ(float z_);

    // Property getting methods
    float GetX() const;
    float GetY() const;
    float GetZ() const;

    // Persistence methods
    virtual void CreateNode(DataNode *node);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, AttributeGroup *rhs);

private:
    float x;
    float y;
    float z;
};

#endif
