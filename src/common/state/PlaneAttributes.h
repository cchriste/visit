#ifndef PLANEATTRIBUTES_H
#define PLANEATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: PlaneAttributes
//
// Purpose:
//    Attributes for a plane
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 13:40:01 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API PlaneAttributes : public AttributeSubject
{
public:
    PlaneAttributes();
    PlaneAttributes(const PlaneAttributes &obj);
    virtual ~PlaneAttributes();

    virtual void operator = (const PlaneAttributes &obj);
    virtual bool operator == (const PlaneAttributes &obj) const;
    virtual bool operator != (const PlaneAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectOrigin();
    void SelectNormal();
    void SelectUpAxis();

    // Property setting methods
    void SetOrigin(const double *origin_);
    void SetNormal(const double *normal_);
    void SetUpAxis(const double *upAxis_);
    void SetHaveRadius(bool haveRadius_);
    void SetRadius(double radius_);

    // Property getting methods
    const double *GetOrigin() const;
          double *GetOrigin();
    const double *GetNormal() const;
          double *GetNormal();
    const double *GetUpAxis() const;
          double *GetUpAxis();
    bool         GetHaveRadius() const;
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
    double normal[3];
    double upAxis[3];
    bool   haveRadius;
    double radius;
};

#endif
