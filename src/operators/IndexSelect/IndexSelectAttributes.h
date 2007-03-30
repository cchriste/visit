#ifndef INDEXSELECTATTRIBUTES_H
#define INDEXSELECTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: IndexSelectAttributes
//
// Purpose:
//    This class contains attributes for the index select operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 14:50:04 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class IndexSelectAttributes : public AttributeSubject
{
public:
    enum Dimension
    {
        OneD,
        TwoD,
        ThreeD
    };
    enum DataType
    {
        AllDomains,
        OneDomain,
        OneGroup
    };

    IndexSelectAttributes();
    IndexSelectAttributes(const IndexSelectAttributes &obj);
    virtual ~IndexSelectAttributes();

    virtual void operator = (const IndexSelectAttributes &obj);
    virtual bool operator == (const IndexSelectAttributes &obj) const;
    virtual bool operator != (const IndexSelectAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetDim(Dimension dim_);
    void SetXMin(int xMin_);
    void SetXMax(int xMax_);
    void SetXIncr(int xIncr_);
    void SetYMin(int yMin_);
    void SetYMax(int yMax_);
    void SetYIncr(int yIncr_);
    void SetZMin(int zMin_);
    void SetZMax(int zMax_);
    void SetZIncr(int zIncr_);
    void SetWhichData(DataType whichData_);
    void SetDomainIndex(int domainIndex_);
    void SetGroupIndex(int groupIndex_);

    // Property getting methods
    Dimension GetDim() const;
    int GetXMin() const;
    int GetXMax() const;
    int GetXIncr() const;
    int GetYMin() const;
    int GetYMax() const;
    int GetYIncr() const;
    int GetZMin() const;
    int GetZMax() const;
    int GetZIncr() const;
    DataType GetWhichData() const;
    int GetDomainIndex() const;
    int GetGroupIndex() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Dimension_ToString(Dimension);
    static bool Dimension_FromString(const std::string &, Dimension &);
protected:
    static std::string Dimension_ToString(int);
public:
    static std::string DataType_ToString(DataType);
    static bool DataType_FromString(const std::string &, DataType &);
protected:
    static std::string DataType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    int dim;
    int xMin;
    int xMax;
    int xIncr;
    int yMin;
    int yMax;
    int yIncr;
    int zMin;
    int zMax;
    int zIncr;
    int whichData;
    int domainIndex;
    int groupIndex;
};

#endif
