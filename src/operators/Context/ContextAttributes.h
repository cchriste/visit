#ifndef CONTEXTATTRIBUTES_H
#define CONTEXTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ContextAttributes
//
// Purpose:
//    This class contains attributes for the context operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Jul 15 14:56:05 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class ContextAttributes : public AttributeSubject
{
public:
    enum Amount
    {
        Some,
        All
    };

    ContextAttributes();
    ContextAttributes(const ContextAttributes &obj);
    virtual ~ContextAttributes();

    virtual void operator = (const ContextAttributes &obj);
    virtual bool operator == (const ContextAttributes &obj) const;
    virtual bool operator != (const ContextAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectContext();

    // Property setting methods
    void SetOffset(double offset_);
    void SetLow(double low_);
    void SetHi(double hi_);
    void SetContext(const std::string &context_);
    void SetCutoff(double cutoff_);
    void SetBelow(double below_);
    void SetAbove(double above_);

    // Property getting methods
    double            GetOffset() const;
    double            GetLow() const;
    double            GetHi() const;
    const std::string &GetContext() const;
          std::string &GetContext();
    double            GetCutoff() const;
    double            GetBelow() const;
    double            GetAbove() const;

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
    double      offset;
    double      low;
    double      hi;
    std::string context;
    double      cutoff;
    double      below;
    double      above;
};

#endif
