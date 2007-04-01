#ifndef STATUSATTRIBUTES_H
#define STATUSATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: StatusAttributes
//
// Purpose:
//    This class contains the status that is displayed in the GUI's status bar.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Mon Nov 1 17:08:56 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class STATE_API StatusAttributes : public AttributeSubject
{
public:
    static const int DEFAULT_DURATION;

    StatusAttributes();
    StatusAttributes(const StatusAttributes &obj);
    virtual ~StatusAttributes();

    virtual void operator = (const StatusAttributes &obj);
    virtual bool operator == (const StatusAttributes &obj) const;
    virtual bool operator != (const StatusAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectSender();
    void SelectStatusMessage();
    void SelectCurrentStageName();

    // Property setting methods
    void SetSender(const std::string &sender_);
    void SetClearStatus(bool clearStatus_);
    void SetStatusMessage(const std::string &statusMessage_);
    void SetPercent(int percent_);
    void SetCurrentStage(int currentStage_);
    void SetCurrentStageName(const std::string &currentStageName_);
    void SetMaxStage(int maxStage_);
    void SetMessageType(int messageType_);
    void SetDuration(int duration_);

    // Property getting methods
    const std::string &GetSender() const;
          std::string &GetSender();
    bool              GetClearStatus() const;
    const std::string &GetStatusMessage() const;
          std::string &GetStatusMessage();
    int               GetPercent() const;
    int               GetCurrentStage() const;
    const std::string &GetCurrentStageName() const;
          std::string &GetCurrentStageName();
    int               GetMaxStage() const;
    int               GetMessageType() const;
    int               GetDuration() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    std::string sender;
    bool        clearStatus;
    std::string statusMessage;
    int         percent;
    int         currentStage;
    std::string currentStageName;
    int         maxStage;
    int         messageType;
    int         duration;
};

#endif
