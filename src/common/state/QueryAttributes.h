#ifndef QUERYATTRIBUTES_H
#define QUERYATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
#include <iostream.h>

// ****************************************************************************
// Class: QueryAttributes
//
// Purpose:
//    This class contains attributes used for query.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 23 11:31:23 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API QueryAttributes : public AttributeSubject
{
public:
    QueryAttributes();
    QueryAttributes(const QueryAttributes &obj);
    virtual ~QueryAttributes();

    virtual void operator = (const QueryAttributes &obj);
    virtual bool operator == (const QueryAttributes &obj) const;
    virtual bool operator != (const QueryAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectName();
    void SelectVariables();
    void SelectResultsMessage();
    void SelectWorldPoint();
    void SelectRayPoint1();
    void SelectRayPoint2();
    void SelectCellPoint();

    // Property setting methods
    void SetName(const std::string &name_);
    void SetVariables(const stringVector &variables_);
    void SetResultsMessage(const std::string &resultsMessage_);
    void SetWorldPoint(const float *worldPoint_);
    void SetDomain(int domain_);
    void SetZone(int zone_);
    void SetRayPoint1(const float *rayPoint1_);
    void SetRayPoint2(const float *rayPoint2_);
    void SetCellPoint(const float *cellPoint_);
    void SetResultsValue(double resultsValue_);

    // Property getting methods
    const std::string  &GetName() const;
          std::string  &GetName();
    const stringVector &GetVariables() const;
          stringVector &GetVariables();
    const std::string  &GetResultsMessage() const;
          std::string  &GetResultsMessage();
    const float        *GetWorldPoint() const;
          float        *GetWorldPoint();
    int                GetDomain() const;
    int                GetZone() const;
    const float        *GetRayPoint1() const;
          float        *GetRayPoint1();
    const float        *GetRayPoint2() const;
          float        *GetRayPoint2();
    const float        *GetCellPoint() const;
          float        *GetCellPoint();
    double             GetResultsValue() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void Reset();
    void PrintSelf(ostream &os);
private:
    std::string  name;
    stringVector variables;
    std::string  resultsMessage;
    float        worldPoint[3];
    int          domain;
    int          zone;
    float        rayPoint1[3];
    float        rayPoint2[3];
    float        cellPoint[3];
    double       resultsValue;
};

#endif
