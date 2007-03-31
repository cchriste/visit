#ifndef BOXEXTENTS_H
#define BOXEXTENTS_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: BoxExtents
//
// Purpose:
//    Attributes for an axis-aligned box
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:23:58 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API BoxExtents : public AttributeSubject
{
public:
    BoxExtents();
    BoxExtents(const BoxExtents &obj);
    virtual ~BoxExtents();

    virtual void operator = (const BoxExtents &obj);
    virtual bool operator == (const BoxExtents &obj) const;
    virtual bool operator != (const BoxExtents &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectExtents();

    // Property setting methods
    void SetExtents(const double *extents_);

    // Property getting methods
    const double *GetExtents() const;
          double *GetExtents();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    double extents[6];
};

#endif
