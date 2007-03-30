#ifndef BOXATTRIBUTES_H
#define BOXATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: BoxAttributes
//
// Purpose:
//    This class contains attributes for the box operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 14:49:41 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class BoxAttributes : public AttributeSubject
{
public:
    enum Amount
    {
        Some,
        All
    };

    BoxAttributes();
    BoxAttributes(const BoxAttributes &obj);
    virtual ~BoxAttributes();

    virtual void operator = (const BoxAttributes &obj);
    virtual bool operator == (const BoxAttributes &obj) const;
    virtual bool operator != (const BoxAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetAmount(Amount amount_);
    void SetMinx(double minx_);
    void SetMaxx(double maxx_);
    void SetMiny(double miny_);
    void SetMaxy(double maxy_);
    void SetMinz(double minz_);
    void SetMaxz(double maxz_);

    // Property getting methods
    Amount GetAmount() const;
    double GetMinx() const;
    double GetMaxx() const;
    double GetMiny() const;
    double GetMaxy() const;
    double GetMinz() const;
    double GetMaxz() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Amount_ToString(Amount);
    static bool Amount_FromString(const std::string &, Amount &);
protected:
    static std::string Amount_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    int    amount;
    double minx;
    double maxx;
    double miny;
    double maxy;
    double minz;
    double maxz;
};

#endif
