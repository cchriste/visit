#ifndef AVTSIMULATIONINFORMATION_H
#define AVTSIMULATIONINFORMATION_H
#include <dbatts_exports.h>
#include <string>
#include <AttributeSubject.h>
class avtSimulationCommandSpecification;

// ****************************************************************************
// Class: avtSimulationInformation
//
// Purpose:
//    Contains information about simulation connections
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Apr 28 11:35:13 PDT 2005
//
// Modifications:
//   
// ****************************************************************************

class DBATTS_API avtSimulationInformation : public AttributeSubject
{
public:
    enum RunMode
    {
        Unknown,
        Running,
        Stopped
    };

    avtSimulationInformation();
    avtSimulationInformation(const avtSimulationInformation &obj);
    virtual ~avtSimulationInformation();

    virtual avtSimulationInformation& operator = (const avtSimulationInformation &obj);
    virtual bool operator == (const avtSimulationInformation &obj) const;
    virtual bool operator != (const avtSimulationInformation &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectHost();
    void SelectSecurityKey();
    void SelectOtherNames();
    void SelectOtherValues();
    void SelectCommands();

    // Property setting methods
    void SetHost(const std::string &host_);
    void SetPort(int port_);
    void SetSecurityKey(const std::string &securityKey_);
    void SetOtherNames(const stringVector &otherNames_);
    void SetOtherValues(const stringVector &otherValues_);
    void SetMode(RunMode mode_);

    // Property getting methods
    const std::string  &GetHost() const;
          std::string  &GetHost();
    int                GetPort() const;
    const std::string  &GetSecurityKey() const;
          std::string  &GetSecurityKey();
    const stringVector &GetOtherNames() const;
          stringVector &GetOtherNames();
    const stringVector &GetOtherValues() const;
          stringVector &GetOtherValues();
    const AttributeGroupVector &GetCommands() const;
          AttributeGroupVector &GetCommands();
    RunMode            GetMode() const;


    // Attributegroup convenience methods
    void AddAvtSimulationCommandSpecification(const avtSimulationCommandSpecification &);
    void ClearAvtSimulationCommandSpecifications();
    void RemoveAvtSimulationCommandSpecification(int i);
    int  GetNumAvtSimulationCommandSpecifications() const;
    avtSimulationCommandSpecification &GetAvtSimulationCommandSpecification(int i);
    const avtSimulationCommandSpecification &GetAvtSimulationCommandSpecification(int i) const;

    avtSimulationCommandSpecification &operator [] (int i);
    const avtSimulationCommandSpecification &operator [] (int i) const;

    // Enum conversion functions
    static std::string RunMode_ToString(RunMode);
    static bool RunMode_FromString(const std::string &, RunMode &);
protected:
    static std::string RunMode_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

protected:
    AttributeGroup *CreateSubAttributeGroup(int index);
private:
    std::string          host;
    int                  port;
    std::string          securityKey;
    stringVector         otherNames;
    stringVector         otherValues;
    AttributeGroupVector commands;
    int                  mode;
};

#endif
