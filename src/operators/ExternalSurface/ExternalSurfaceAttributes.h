#ifndef EXTERNALSURFACEATTRIBUTES_H
#define EXTERNALSURFACEATTRIBUTES_H
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ExternalSurfaceAttributes
//
// Purpose:
//    This class contains attributes for the external surface operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Aug 16 14:32:12 PST 2005
//
// Modifications:
//   
// ****************************************************************************

class ExternalSurfaceAttributes : public AttributeSubject
{
public:
    ExternalSurfaceAttributes();
    ExternalSurfaceAttributes(const ExternalSurfaceAttributes &obj);
    virtual ~ExternalSurfaceAttributes();

    virtual ExternalSurfaceAttributes& operator = (const ExternalSurfaceAttributes &obj);
    virtual bool operator == (const ExternalSurfaceAttributes &obj) const;
    virtual bool operator != (const ExternalSurfaceAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetDummy(bool dummy_);

    // Property getting methods
    bool GetDummy() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool dummy;
};

#endif
