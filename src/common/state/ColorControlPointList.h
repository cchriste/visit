#ifndef COLORCONTROLPOINTLIST_H
#define COLORCONTROLPOINTLIST_H
#include <state_exports.h>
#include <AttributeSubject.h>
class ColorControlPoint;

// ****************************************************************************
// Class: ColorControlPointList
//
// Purpose:
//    This class contains a list of ColorControlPoint objects.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 23 11:29:47 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ColorControlPointList : public AttributeSubject
{
public:
    ColorControlPointList();
    ColorControlPointList(const ColorControlPointList &obj);
    virtual ~ColorControlPointList();

    virtual void operator = (const ColorControlPointList &obj);
    virtual bool operator == (const ColorControlPointList &obj) const;
    virtual bool operator != (const ColorControlPointList &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectControlPoints();

    // Property setting methods
    void SetControlPoints(const AttributeGroupVector &controlPoints_);
    void SetSmoothingFlag(bool smoothingFlag_);
    void SetEqualSpacingFlag(bool equalSpacingFlag_);
    void SetDiscreteFlag(bool discreteFlag_);
    void SetExternalFlag(bool externalFlag_);

    // Property getting methods
    const AttributeGroupVector &GetControlPoints() const;
          AttributeGroupVector &GetControlPoints();
    bool GetSmoothingFlag() const;
    bool GetEqualSpacingFlag() const;
    bool GetDiscreteFlag() const;
    bool GetExternalFlag() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddColorControlPoint(const ColorControlPoint &);
    void ClearColorControlPoints();
    void RemoveColorControlPoint(int i);
    int  GetNumColorControlPoints() const;
    ColorControlPoint &GetColorControlPoint(int i);
    const ColorControlPoint &GetColorControlPoint(int i) const;

    ColorControlPoint &operator [] (int i);
    const ColorControlPoint &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void GetColors(unsigned char *rgb, int ncolors) const;
protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    AttributeGroupVector controlPoints;
    bool                 smoothingFlag;
    bool                 equalSpacingFlag;
    bool                 discreteFlag;
    bool                 externalFlag;
};

#endif
