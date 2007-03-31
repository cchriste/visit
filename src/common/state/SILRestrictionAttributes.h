#ifndef SILRESTRICTIONATTRIBUTES_H
#define SILRESTRICTIONATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>
#include <SILAttributes.h>

// ****************************************************************************
// Class: SILRestrictionAttributes
//
// Purpose:
//    The class contains attributes for SIL restrictions.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:24:30 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SILRestrictionAttributes : public AttributeSubject
{
public:
    SILRestrictionAttributes();
    SILRestrictionAttributes(const SILRestrictionAttributes &obj);
    virtual ~SILRestrictionAttributes();

    virtual void operator = (const SILRestrictionAttributes &obj);
    virtual bool operator == (const SILRestrictionAttributes &obj) const;
    virtual bool operator != (const SILRestrictionAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectUseSet();
    void SelectSilAtts();

    // Property setting methods
    void SetUseSet(const unsignedCharVector &useSet_);
    void SetTopSet(int topSet_);
    void SetSilAtts(const SILAttributes &silAtts_);

    // Property getting methods
    const unsignedCharVector &GetUseSet() const;
          unsignedCharVector &GetUseSet();
    int                      GetTopSet() const;
    const SILAttributes      &GetSilAtts() const;
          SILAttributes      &GetSilAtts();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    unsignedCharVector useSet;
    int                topSet;
    SILAttributes      silAtts;
};

#endif
