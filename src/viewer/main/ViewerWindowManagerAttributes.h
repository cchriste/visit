#ifndef VIEWERWINDOWMANAGERATTRIBUTES_H
#define VIEWERWINDOWMANAGERATTRIBUTES_H
#include <AttributeSubject.h>
class ActionGroupDescription;
#include "ViewerRPC.h"

// ****************************************************************************
// Class: ViewerWindowManagerAttributes
//
// Purpose:
//    This class contains the attributes that dictate where viewer windows are positioned, etc.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Mar 16 13:10:38 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class ViewerWindowManagerAttributes : public AttributeSubject
{
public:
    ViewerWindowManagerAttributes();
    ViewerWindowManagerAttributes(const ViewerWindowManagerAttributes &obj);
    virtual ~ViewerWindowManagerAttributes();

    virtual void operator = (const ViewerWindowManagerAttributes &obj);
    virtual bool operator == (const ViewerWindowManagerAttributes &obj) const;
    virtual bool operator != (const ViewerWindowManagerAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectActionConfigurations();

    // Property setting methods
    void SetActionConfigurations(const AttributeGroupVector &actionConfigurations_);
    void SetToolbarsVisible(bool toolbarsVisible_);
    void SetLargeIcons(bool largeIcons_);

    // Property getting methods
    const AttributeGroupVector &GetActionConfigurations() const;
          AttributeGroupVector &GetActionConfigurations();
    bool GetToolbarsVisible() const;
    bool GetLargeIcons() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Attributegroup convenience methods
    void AddActionGroupDescription(const ActionGroupDescription &);
    void ClearActionGroupDescriptions();
    void RemoveActionGroupDescription(int i);
    int  GetNumActionGroupDescriptions() const;
    ActionGroupDescription &GetActionGroupDescription(int i);
    const ActionGroupDescription &GetActionGroupDescription(int i) const;

    ActionGroupDescription &operator [] (int i);
    const ActionGroupDescription &operator [] (int i) const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual void ProcessOldVersions(DataNode *parentNode, const char *configVersion);
    void RemoveActionFromNode(DataNode *, const char *, ViewerRPC::ViewerRPCType);
    void AddAction(DataNode *, const char *, ViewerRPC::ViewerRPCType);
    void AddActionGroup(DataNode *, ActionGroupDescription &);
protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    AttributeGroupVector actionConfigurations;
    bool                 toolbarsVisible;
    bool                 largeIcons;
};

#endif
