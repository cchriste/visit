#ifndef SYNCATTRIBUTES_H
#define SYNCATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: SyncAttributes
//
// Purpose:
//    This class contains an integer that can be used to synchronize the viewer and its clients.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 23 11:31:57 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SyncAttributes : public AttributeSubject
{
public:
    SyncAttributes();
    SyncAttributes(const SyncAttributes &obj);
    virtual ~SyncAttributes();

    virtual void operator = (const SyncAttributes &obj);
    virtual bool operator == (const SyncAttributes &obj) const;
    virtual bool operator != (const SyncAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetSyncTag(int syncTag_);

    // Property getting methods
    int GetSyncTag() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    int syncTag;
};

#endif
