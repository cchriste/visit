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
// Creation:   Tue Jun 24 11:35:39 PDT 2003
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

    // Property selection methods
    virtual void SelectAll();
    void SelectWindowX();
    void SelectWindowY();
    void SelectWindowWidth();
    void SelectWindowHeight();
    void SelectActionConfigurations();

    // Property setting methods
    void SetNWindows(int nWindows_);
    void SetWindowX(const intVector &windowX_);
    void SetWindowY(const intVector &windowY_);
    void SetWindowWidth(const intVector &windowWidth_);
    void SetWindowHeight(const intVector &windowHeight_);
    void SetActionConfigurations(const AttributeGroupVector &actionConfigurations_);

    // Property getting methods
    int             GetNWindows() const;
    const intVector &GetWindowX() const;
          intVector &GetWindowX();
    const intVector &GetWindowY() const;
          intVector &GetWindowY();
    const intVector &GetWindowWidth() const;
          intVector &GetWindowWidth();
    const intVector &GetWindowHeight() const;
          intVector &GetWindowHeight();
    const AttributeGroupVector &GetActionConfigurations() const;
          AttributeGroupVector &GetActionConfigurations();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
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
    int                  nWindows;
    intVector            windowX;
    intVector            windowY;
    intVector            windowWidth;
    intVector            windowHeight;
    AttributeGroupVector actionConfigurations;
};

#endif
