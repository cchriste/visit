#ifndef MINMAXINFO_H
#define MINMAXINFO_H
#include <query_exports.h>
#include <string>
#include <AttributeSubject.h>
#include <avtMatrix.h>

// ****************************************************************************
// Class: MinMaxInfo
//
// Purpose:
//    Results storage for MinMaxQuery.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Wed Jul 28 08:58:20 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

class QUERY_API MinMaxInfo : public AttributeSubject
{
public:
    MinMaxInfo();
    MinMaxInfo(const MinMaxInfo &obj);
    virtual ~MinMaxInfo();

    virtual void operator = (const MinMaxInfo &obj);
    virtual bool operator == (const MinMaxInfo &obj) const;
    virtual bool operator != (const MinMaxInfo &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectCoord();
    void SelectType();
    void SelectMatName();

    // Property setting methods
    void SetElementNum(int elementNum_);
    void SetDomain(int domain_);
    void SetValue(float value_);
    void SetCoord(const float *coord_);
    void SetType(const std::string &type_);
    void SetMatName(const std::string &matName_);

    // Property getting methods
    int               GetElementNum() const;
    int               GetDomain() const;
    float             GetValue() const;
    const float       *GetCoord() const;
          float       *GetCoord();
    const std::string &GetType() const;
          std::string &GetType();
    const std::string &GetMatName() const;
          std::string &GetMatName();


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void Initialize(const float, const std::string &);
    void TransformCoord(const avtMatrix *trans);
    bool  EquivalentForOutput (const MinMaxInfo &obj) const;
private:
    int         elementNum;
    int         domain;
    float       value;
    float       coord[3];
    std::string type;
    std::string matName;
};

#endif
