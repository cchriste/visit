#ifndef GAUSSIANCONTROLPOINTLIST_H
#define GAUSSIANCONTROLPOINTLIST_H
#include <state_exports.h>
#include <AttributeSubject.h>
class GaussianControlPoint;

// ****************************************************************************
// Class: GaussianControlPointList
//
// Purpose:
//    This class contains a list of GaussianControlPoint objects.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 23 11:30:12 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API GaussianControlPointList : public AttributeSubject
{
public:
    GaussianControlPointList();
    GaussianControlPointList(const GaussianControlPointList &obj);
    virtual ~GaussianControlPointList();

    virtual void operator = (const GaussianControlPointList &obj);
    virtual bool operator == (const GaussianControlPointList &obj) const;
    virtual bool operator != (const GaussianControlPointList &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectControlPoints();

    // Property setting methods
    void SetControlPoints(const AttributeGroupVector &controlPoints_);

    // Property getting methods
    const AttributeGroupVector &GetControlPoints() const;
          AttributeGroupVector &GetControlPoints();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddGaussianControlPoint(const GaussianControlPoint &);
    void ClearGaussianControlPoints();
    void RemoveGaussianControlPoint(int i);
    int  GetNumGaussianControlPoints() const;
    GaussianControlPoint &GetGaussianControlPoint(int i);
    const GaussianControlPoint &GetGaussianControlPoint(int i) const;

    GaussianControlPoint &operator [] (int i);
    const GaussianControlPoint &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    AttributeGroupVector controlPoints;
};

#endif
