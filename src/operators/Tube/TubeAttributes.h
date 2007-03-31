#ifndef TUBEATTRIBUTES_H
#define TUBEATTRIBUTES_H
#include <AttributeSubject.h>

// ****************************************************************************
// Class: TubeAttributes
//
// Purpose:
//    This class contains attributes for the tube operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:50:27 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class TubeAttributes : public AttributeSubject
{
public:
    TubeAttributes();
    TubeAttributes(const TubeAttributes &obj);
    virtual ~TubeAttributes();

    virtual void operator = (const TubeAttributes &obj);
    virtual bool operator == (const TubeAttributes &obj) const;
    virtual bool operator != (const TubeAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetWidth(float width_);
    void SetFineness(int fineness_);
    void SetCapping(bool capping_);

    // Property getting methods
    float GetWidth() const;
    int   GetFineness() const;
    bool  GetCapping() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    float width;
    int   fineness;
    bool  capping;
};

#endif
