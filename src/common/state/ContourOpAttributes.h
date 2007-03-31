#ifndef CONTOUROPATTRIBUTES_H
#define CONTOUROPATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ContourOpAttributes
//
// Purpose:
//    This class contains the operator attributes for the contour operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Jul 15 13:47:30 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ContourOpAttributes : public AttributeSubject
{
public:
    enum ContourMethod
    {
        Level,
        Value,
        Percent
    };
    enum ContourScaling
    {
        Linear,
        Log
    };

    ContourOpAttributes();
    ContourOpAttributes(const ContourOpAttributes &obj);
    virtual ~ContourOpAttributes();

    virtual void operator = (const ContourOpAttributes &obj);
    virtual bool operator == (const ContourOpAttributes &obj) const;
    virtual bool operator != (const ContourOpAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectContourValue();
    void SelectContourPercent();
    void SelectVariable();

    // Property setting methods
    void SetContourNLevels(int contourNLevels_);
    void SetContourValue(const doubleVector &contourValue_);
    void SetContourPercent(const doubleVector &contourPercent_);
    void SetContourMethod(ContourMethod contourMethod_);
    void SetMinFlag(bool minFlag_);
    void SetMaxFlag(bool maxFlag_);
    void SetMin(double min_);
    void SetMax(double max_);
    void SetScaling(ContourScaling scaling_);
    void SetVariable(const std::string &variable_);

    // Property getting methods
    int                GetContourNLevels() const;
    const doubleVector &GetContourValue() const;
          doubleVector &GetContourValue();
    const doubleVector &GetContourPercent() const;
          doubleVector &GetContourPercent();
    ContourMethod      GetContourMethod() const;
    bool               GetMinFlag() const;
    bool               GetMaxFlag() const;
    double             GetMin() const;
    double             GetMax() const;
    ContourScaling     GetScaling() const;
    const std::string  &GetVariable() const;
          std::string  &GetVariable();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string ContourMethod_ToString(ContourMethod);
    static bool ContourMethod_FromString(const std::string &, ContourMethod &);
protected:
    static std::string ContourMethod_ToString(int);
public:
    static std::string ContourScaling_ToString(ContourScaling);
    static bool ContourScaling_FromString(const std::string &, ContourScaling &);
protected:
    static std::string ContourScaling_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    int          contourNLevels;
    doubleVector contourValue;
    doubleVector contourPercent;
    int          contourMethod;
    bool         minFlag;
    bool         maxFlag;
    double       min;
    double       max;
    int          scaling;
    std::string  variable;
};

#endif
