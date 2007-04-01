#ifndef GLOBALLINEOUTATTRIBUTES_H
#define GLOBALLINEOUTATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: GlobalLineoutAttributes
//
// Purpose:
//    This file contains global attributes controlling Lineouts.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Fri Nov 19 10:43:37 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

class STATE_API GlobalLineoutAttributes : public AttributeSubject
{
public:
    GlobalLineoutAttributes();
    GlobalLineoutAttributes(const GlobalLineoutAttributes &obj);
    virtual ~GlobalLineoutAttributes();

    virtual GlobalLineoutAttributes& operator = (const GlobalLineoutAttributes &obj);
    virtual bool operator == (const GlobalLineoutAttributes &obj) const;
    virtual bool operator != (const GlobalLineoutAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetDynamic(bool Dynamic_);
    void SetCreateWindow(bool createWindow_);
    void SetWindowId(int windowId_);
    void SetSamplingOn(bool samplingOn_);
    void SetNumSamples(int numSamples_);
    void SetCreateReflineLabels(bool createReflineLabels_);

    // Property getting methods
    bool GetDynamic() const;
    bool GetCreateWindow() const;
    int  GetWindowId() const;
    bool GetSamplingOn() const;
    int  GetNumSamples() const;
    bool GetCreateReflineLabels() const;


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool Dynamic;
    bool createWindow;
    int  windowId;
    bool samplingOn;
    int  numSamples;
    bool createReflineLabels;
};

#endif
