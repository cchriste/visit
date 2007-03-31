#ifndef CONNCOMPREDUCEATTRIBUTES_H
#define CONNCOMPREDUCEATTRIBUTES_H
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ConnCompReduceAttributes
//
// Purpose:
//    This class contains attributes for the reduce connected components operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:49:29 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class ConnCompReduceAttributes : public AttributeSubject
{
public:
    ConnCompReduceAttributes();
    ConnCompReduceAttributes(const ConnCompReduceAttributes &obj);
    virtual ~ConnCompReduceAttributes();

    virtual void operator = (const ConnCompReduceAttributes &obj);
    virtual bool operator == (const ConnCompReduceAttributes &obj) const;
    virtual bool operator != (const ConnCompReduceAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetTarget(double target_);

    // Property getting methods
    double GetTarget() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    double target;
};

#endif
