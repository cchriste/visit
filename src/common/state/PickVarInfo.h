#ifndef PICKVARINFO_H
#define PICKVARINFO_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
#include "snprintf.h"

// ****************************************************************************
// Class: PickVarInfo
//
// Purpose:
//    This class contains PickVarInfo.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Jul 15 13:47:42 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API PickVarInfo : public AttributeSubject
{
public:
    enum Centering
    {
        Nodal,
        Zonal,
        None
    };

    PickVarInfo();
    PickVarInfo(const PickVarInfo &obj);
    virtual ~PickVarInfo();

    virtual void operator = (const PickVarInfo &obj);
    virtual bool operator == (const PickVarInfo &obj) const;
    virtual bool operator != (const PickVarInfo &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVariableName();
    void SelectNames();
    void SelectValues();
    void SelectMixNames();
    void SelectMixValues();
    void SelectMiscMessage();
    void SelectNumMatsPerZone();

    // Property setting methods
    void SetVariableName(const std::string &variableName_);
    void SetNames(const stringVector &names_);
    void SetValues(const doubleVector &values_);
    void SetMixNames(const stringVector &mixNames_);
    void SetMixValues(const doubleVector &mixValues_);
    void SetMixVar(bool mixVar_);
    void SetCentering(Centering centering_);
    void SetMiscMessage(const std::string &miscMessage_);
    void SetVarIsMaterial(bool varIsMaterial_);
    void SetNumMatsPerZone(const intVector &numMatsPerZone_);

    // Property getting methods
    const std::string  &GetVariableName() const;
          std::string  &GetVariableName();
    const stringVector &GetNames() const;
          stringVector &GetNames();
    const doubleVector &GetValues() const;
          doubleVector &GetValues();
    const stringVector &GetMixNames() const;
          stringVector &GetMixNames();
    const doubleVector &GetMixValues() const;
          doubleVector &GetMixValues();
    bool               GetMixVar() const;
    Centering          GetCentering() const;
    const std::string  &GetMiscMessage() const;
          std::string  &GetMiscMessage();
    bool               GetVarIsMaterial() const;
    const intVector    &GetNumMatsPerZone() const;
          intVector    &GetNumMatsPerZone();

    // Enum conversion functions
    static std::string Centering_ToString(Centering);
    static bool Centering_FromString(const std::string &, Centering &);
protected:
    static std::string Centering_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void PrintSelf(ostream &os);
    void CreateOutputString(std::string &os);
    void CreateOutputStrings(std::vector<std::string> &os);
private:
    std::string  variableName;
    stringVector names;
    doubleVector values;
    stringVector mixNames;
    doubleVector mixValues;
    bool         mixVar;
    int          centering;
    std::string  miscMessage;
    bool         varIsMaterial;
    intVector    numMatsPerZone;
};

#endif
