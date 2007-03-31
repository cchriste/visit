#ifndef VIEWCURVEATTRIBUTES_H
#define VIEWCURVEATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ViewCurveAttributes
//
// Purpose:
//    This class contains the curve view attributes.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:24:38 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ViewCurveAttributes : public AttributeSubject
{
public:
    ViewCurveAttributes();
    ViewCurveAttributes(const ViewCurveAttributes &obj);
    virtual ~ViewCurveAttributes();

    virtual void operator = (const ViewCurveAttributes &obj);
    virtual bool operator == (const ViewCurveAttributes &obj) const;
    virtual bool operator != (const ViewCurveAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectDomainCoords();
    void SelectRangeCoords();
    void SelectViewportCoords();

    // Property setting methods
    void SetDomainCoords(const double *domainCoords_);
    void SetRangeCoords(const double *rangeCoords_);
    void SetViewportCoords(const double *viewportCoords_);

    // Property getting methods
    const double *GetDomainCoords() const;
          double *GetDomainCoords();
    const double *GetRangeCoords() const;
          double *GetRangeCoords();
    const double *GetViewportCoords() const;
          double *GetViewportCoords();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    double domainCoords[2];
    double rangeCoords[2];
    double viewportCoords[4];
};

#endif
