#ifndef GAUSSIANCONTROLPOINT_H
#define GAUSSIANCONTROLPOINT_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: GaussianControlPoint
//
// Purpose:
//    This class contains the information for a gaussian in the opacity bar.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 13:39:50 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API GaussianControlPoint : public AttributeSubject
{
public:
    GaussianControlPoint();
    GaussianControlPoint(const GaussianControlPoint &obj);
    virtual ~GaussianControlPoint();

    virtual void operator = (const GaussianControlPoint &obj);
    virtual bool operator == (const GaussianControlPoint &obj) const;
    virtual bool operator != (const GaussianControlPoint &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetX(float x_);
    void SetHeight(float height_);
    void SetWidth(float width_);
    void SetXBias(float xBias_);
    void SetYBias(float yBias_);

    // Property getting methods
    float GetX() const;
    float GetHeight() const;
    float GetWidth() const;
    float GetXBias() const;
    float GetYBias() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    float x;
    float height;
    float width;
    float xBias;
    float yBias;
};

#endif
