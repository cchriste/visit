#ifndef INVERSEGHOSTZONEATTRIBUTES_H
#define INVERSEGHOSTZONEATTRIBUTES_H
#include <AttributeSubject.h>

// ****************************************************************************
// Class: InverseGhostZoneAttributes
//
// Purpose:
//    This class contains attributes for the inverse ghost zone operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:49:47 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class InverseGhostZoneAttributes : public AttributeSubject
{
public:
    InverseGhostZoneAttributes();
    InverseGhostZoneAttributes(const InverseGhostZoneAttributes &obj);
    virtual ~InverseGhostZoneAttributes();

    virtual void operator = (const InverseGhostZoneAttributes &obj);
    virtual bool operator == (const InverseGhostZoneAttributes &obj) const;
    virtual bool operator != (const InverseGhostZoneAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetConstantData(bool constantData_);

    // Property getting methods
    bool GetConstantData() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool constantData;
};

#endif
