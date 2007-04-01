#ifndef INTERACTORATTRIBUTES_H
#define INTERACTORATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: InteractorAttributes
//
// Purpose:
//    This class contains attributes associated with the main window.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Nov 11 14:04:09 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class STATE_API InteractorAttributes : public AttributeSubject
{
public:
    enum NavigationMode
    {
        Trackball,
        Flythrough
    };

    InteractorAttributes();
    InteractorAttributes(const InteractorAttributes &obj);
    virtual ~InteractorAttributes();

    virtual InteractorAttributes& operator = (const InteractorAttributes &obj);
    virtual bool operator == (const InteractorAttributes &obj) const;
    virtual bool operator != (const InteractorAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetShowGuidelines(bool showGuidelines_);
    void SetClampSquare(bool clampSquare_);
    void SetNavigationMode(NavigationMode navigationMode_);

    // Property getting methods
    bool GetShowGuidelines() const;
    bool GetClampSquare() const;
    NavigationMode GetNavigationMode() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string NavigationMode_ToString(NavigationMode);
    static bool NavigationMode_FromString(const std::string &, NavigationMode &);
protected:
    static std::string NavigationMode_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool showGuidelines;
    bool clampSquare;
    int  navigationMode;
};

#endif
