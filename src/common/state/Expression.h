#ifndef EXPRESSION_H
#define EXPRESSION_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: Expression
//
// Purpose:
//    This class contains an expression.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Jul 24 11:03:22 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API Expression : public AttributeSubject
{
public:
    enum ExprType
    {
        Unknown,
        ScalarMeshVar,
        VectorMeshVar,
        TensorMeshVar,
        SymmetricTensorMeshVar,
        Mesh,
        Material,
        Species
    };

    Expression();
    Expression(const Expression &obj);
    virtual ~Expression();

    virtual void operator = (const Expression &obj);
    virtual bool operator == (const Expression &obj) const;
    virtual bool operator != (const Expression &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectName();
    void SelectDefinition();

    // Property setting methods
    void SetName(const std::string &name_);
    void SetDefinition(const std::string &definition_);
    void SetHidden(bool hidden_);
    void SetType(ExprType type_);

    // Property getting methods
    const std::string &GetName() const;
          std::string &GetName();
    const std::string &GetDefinition() const;
          std::string &GetDefinition();
    bool              GetHidden() const;
    ExprType          GetType() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string ExprType_ToString(ExprType);
    static bool ExprType_FromString(const std::string &, ExprType &);
protected:
    static std::string ExprType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    static const int GetNumTypes();
    static const char * GetTypeString(const Expression::ExprType t);
    static const Expression::ExprType GetTypeId(const std::string s);
private:
    std::string name;
    std::string definition;
    bool        hidden;
    int         type;
};

#endif
