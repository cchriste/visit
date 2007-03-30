#ifndef SPHEREATTRIBUTES_H
#define SPHEREATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: SphereAttributes
//
// Purpose:
//    Attributes for a sphere
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 13:40:12 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SphereAttributes : public AttributeSubject
{
public:
    SphereAttributes();
    SphereAttributes(const SphereAttributes &obj);
    virtual ~SphereAttributes();

    virtual void operator = (const SphereAttributes &obj);
    virtual bool operator == (const SphereAttributes &obj) const;
    virtual bool operator != (const SphereAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectOrigin();

    // Property setting methods
    void SetOrigin(const double *origin_);
    void SetRadius(double radius_);

    // Property getting methods
    const double *GetOrigin() const;
          double *GetOrigin();
    double       GetRadius() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    double origin[3];
    double radius;
};

#endif
