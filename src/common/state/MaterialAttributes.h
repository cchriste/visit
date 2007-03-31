#ifndef MATERIALATTRIBUTES_H
#define MATERIALATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: MaterialAttributes
//
// Purpose:
//    Attributes to control material interface reconstruction
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Jul 29 12:56:51 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API MaterialAttributes : public AttributeSubject
{
public:
    MaterialAttributes();
    MaterialAttributes(const MaterialAttributes &obj);
    virtual ~MaterialAttributes();

    virtual void operator = (const MaterialAttributes &obj);
    virtual bool operator == (const MaterialAttributes &obj) const;
    virtual bool operator != (const MaterialAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetSmoothing(bool smoothing_);
    void SetForceMIR(bool forceMIR_);
    void SetCleanZonesOnly(bool cleanZonesOnly_);
    void SetNeedValidConnectivity(bool needValidConnectivity_);

    // Property getting methods
    bool GetSmoothing() const;
    bool GetForceMIR() const;
    bool GetCleanZonesOnly() const;
    bool GetNeedValidConnectivity() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool smoothing;
    bool forceMIR;
    bool cleanZonesOnly;
    bool needValidConnectivity;
};

#endif
