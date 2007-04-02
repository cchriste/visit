#ifndef DBPLUGININFOATTRIBUTES_H
#define DBPLUGININFOATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
class DBOptionsAttributes;

// ****************************************************************************
// Class: DBPluginInfoAttributes
//
// Purpose:
//    This class contains the attributes for all the database plugins.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu May 26 14:25:37 PST 2005
//
// Modifications:
//   
// ****************************************************************************

class STATE_API DBPluginInfoAttributes : public AttributeSubject
{
public:
    DBPluginInfoAttributes();
    DBPluginInfoAttributes(const DBPluginInfoAttributes &obj);
    virtual ~DBPluginInfoAttributes();

    virtual DBPluginInfoAttributes& operator = (const DBPluginInfoAttributes &obj);
    virtual bool operator == (const DBPluginInfoAttributes &obj) const;
    virtual bool operator != (const DBPluginInfoAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectTypes();
    void SelectHasWriter();
    void SelectDbOptions();
    void SelectTypesFullNames();

    // Property setting methods
    void SetTypes(const stringVector &types_);
    void SetHasWriter(const intVector &hasWriter_);
    void SetTypesFullNames(const stringVector &typesFullNames_);

    // Property getting methods
    const stringVector &GetTypes() const;
          stringVector &GetTypes();
    const intVector    &GetHasWriter() const;
          intVector    &GetHasWriter();
    const AttributeGroupVector &GetDbOptions() const;
          AttributeGroupVector &GetDbOptions();
    const stringVector &GetTypesFullNames() const;
          stringVector &GetTypesFullNames();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddDBOptionsAttributes(const DBOptionsAttributes &);
    void ClearDBOptionsAttributess();
    void RemoveDBOptionsAttributes(int i);
    int  GetNumDBOptionsAttributess() const;
    DBOptionsAttributes &GetDBOptionsAttributes(int i);
    const DBOptionsAttributes &GetDBOptionsAttributes(int i) const;

    DBOptionsAttributes &operator [] (int i);
    const DBOptionsAttributes &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    stringVector         types;
    intVector            hasWriter;
    AttributeGroupVector dbOptions;
    stringVector         typesFullNames;
};

#endif
