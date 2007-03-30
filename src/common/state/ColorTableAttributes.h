#ifndef COLORTABLEATTRIBUTES_H
#define COLORTABLEATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
class ColorControlPointList;

// ****************************************************************************
// Class: ColorTableAttributes
//
// Purpose:
//    This class contains the list of colortables.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 13:39:45 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API ColorTableAttributes : public AttributeSubject
{
public:
    ColorTableAttributes();
    ColorTableAttributes(const ColorTableAttributes &obj);
    virtual ~ColorTableAttributes();

    virtual void operator = (const ColorTableAttributes &obj);
    virtual bool operator == (const ColorTableAttributes &obj) const;
    virtual bool operator != (const ColorTableAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectNames();
    void SelectColorTables();
    void SelectActiveContinuous();
    void SelectActiveDiscrete();

    // Property setting methods
    void SetNames(const stringVector &names_);
    void SetColorTables(const AttributeGroupVector &colorTables_);
    void SetActiveContinuous(const std::string &activeContinuous_);
    void SetActiveDiscrete(const std::string &activeDiscrete_);

    // Property getting methods
    const stringVector &GetNames() const;
          stringVector &GetNames();
    const AttributeGroupVector &GetColorTables() const;
          AttributeGroupVector &GetColorTables();
    const std::string  &GetActiveContinuous() const;
          std::string  &GetActiveContinuous();
    const std::string  &GetActiveDiscrete() const;
          std::string  &GetActiveDiscrete();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddColorControlPointList(const ColorControlPointList &);
    void ClearColorControlPointLists();
    void RemoveColorControlPointList(int i);
    int  GetNumColorControlPointLists() const;
    ColorControlPointList &GetColorControlPointList(int i);
    const ColorControlPointList &GetColorControlPointList(int i) const;

    ColorControlPointList &operator [] (int i);
    const ColorControlPointList &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    int GetColorTableIndex(const std::string &name) const;
    const ColorControlPointList *GetColorControlPoints(int index) const;
    const ColorControlPointList *GetColorControlPoints(const std::string &name) const;
    void AddColorTable(const std::string &name, const ColorControlPointList &cpts);
    void RemoveColorTable(const std::string &name);
    void RemoveColorTable(int index);
    int GetNumColorTables() const;
protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    stringVector         names;
    AttributeGroupVector colorTables;
    std::string          activeContinuous;
    std::string          activeDiscrete;
};

#endif
