#ifndef LINEOUTATTRIBUTES_H
#define LINEOUTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>
#include <visitstream.h>

// ****************************************************************************
// Class: LineoutAttributes
//
// Purpose:
//    Attributes for the Lineout operator 
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 28 14:34:34 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class LineoutAttributes : public AttributeSubject
{
public:
    LineoutAttributes();
    LineoutAttributes(const LineoutAttributes &obj);
    virtual ~LineoutAttributes();

    virtual void operator = (const LineoutAttributes &obj);
    virtual bool operator == (const LineoutAttributes &obj) const;
    virtual bool operator != (const LineoutAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPoint1();
    void SelectPoint2();
    void SelectDesignator();

    // Property setting methods
    void SetPoint1(const double *point1_);
    void SetPoint2(const double *point2_);
    void SetSamplingOn(bool samplingOn_);
    void SetNumberOfSamplePoints(int numberOfSamplePoints_);
    void SetInteractive(bool interactive_);
    void SetReflineLabels(bool reflineLabels_);
    void SetDesignator(const std::string &designator_);

    // Property getting methods
    const double      *GetPoint1() const;
          double      *GetPoint1();
    const double      *GetPoint2() const;
          double      *GetPoint2();
    bool              GetSamplingOn() const;
    int               GetNumberOfSamplePoints() const;
    bool              GetInteractive() const;
    bool              GetReflineLabels() const;
    const std::string &GetDesignator() const;
          std::string &GetDesignator();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void Print(ostream &, bool) const;
private:
    double      point1[3];
    double      point2[3];
    bool        samplingOn;
    int         numberOfSamplePoints;
    bool        interactive;
    bool        reflineLabels;
    std::string designator;
};

#endif
