#ifndef DATABASECORRELATION_H
#define DATABASECORRELATION_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>
#include <visitstream.h>

// ****************************************************************************
// Class: DatabaseCorrelation
//
// Purpose:
//    This class encapsulates a database correlation, which is a mapping of one or more databases to a set of indices that go from 0 to N.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Feb 3 19:04:13 PST 2004
//
// Modifications:
//   
// ****************************************************************************

class STATE_API DatabaseCorrelation : public AttributeSubject
{
public:
    enum CorrelationMethod
    {
        IndexForIndexCorrelation,
        StretchedIndexCorrelation,
        TimeCorrelation,
        CycleCorrelation,
        UserDefinedCorrelation
    };

    DatabaseCorrelation();
    DatabaseCorrelation(const DatabaseCorrelation &obj);
    virtual ~DatabaseCorrelation();

    virtual void operator = (const DatabaseCorrelation &obj);
    virtual bool operator == (const DatabaseCorrelation &obj) const;
    virtual bool operator != (const DatabaseCorrelation &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectName();
    void SelectDatabaseNames();
    void SelectDatabaseNStates();
    void SelectDatabaseTimes();
    void SelectDatabaseCycles();
    void SelectIndices();
    void SelectCondensedTimes();
    void SelectCondensedCycles();

    // Property setting methods
    void SetName(const std::string &name_);
    void SetNumStates(int numStates_);
    void SetMethod(CorrelationMethod method_);
    void SetDatabaseNames(const stringVector &databaseNames_);
    void SetDatabaseNStates(const intVector &databaseNStates_);
    void SetDatabaseTimes(const doubleVector &databaseTimes_);
    void SetDatabaseCycles(const intVector &databaseCycles_);
    void SetIndices(const intVector &indices_);
    void SetCondensedTimes(const doubleVector &condensedTimes_);
    void SetCondensedCycles(const intVector &condensedCycles_);

    // Property getting methods
    const std::string  &GetName() const;
          std::string  &GetName();
    int                GetNumStates() const;
    CorrelationMethod  GetMethod() const;
    const stringVector &GetDatabaseNames() const;
          stringVector &GetDatabaseNames();
    const intVector    &GetDatabaseNStates() const;
          intVector    &GetDatabaseNStates();
    const doubleVector &GetDatabaseTimes() const;
          doubleVector &GetDatabaseTimes();
    const intVector    &GetDatabaseCycles() const;
          intVector    &GetDatabaseCycles();
    const intVector    &GetIndices() const;
          intVector    &GetIndices();
    const doubleVector &GetCondensedTimes() const;
          doubleVector &GetCondensedTimes();
    const intVector    &GetCondensedCycles() const;
          intVector    &GetCondensedCycles();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string CorrelationMethod_ToString(CorrelationMethod);
    static bool CorrelationMethod_FromString(const std::string &, CorrelationMethod &);
protected:
    static std::string CorrelationMethod_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool GetCorrelatedTimeStates(int state, intVector &states) const;
    void AddDatabase(const std::string &database, int nStates, const doubleVector &times, const intVector &cycles);
    bool UsesDatabase(const std::string &database) const;
    int GetNumDatabases() const;
    int GetCorrelatedTimeState(const std::string &db, int state) const;
    int GetInverseCorrelatedTimeState(const std::string &db, int state) const;
    int GetCondensedCycleForState(int state) const;
    double GetCondensedTimeForState(int state) const;
private:
    std::string  name;
    int          numStates;
    int          method;
    stringVector databaseNames;
    intVector    databaseNStates;
    doubleVector databaseTimes;
    intVector    databaseCycles;
    intVector    indices;
    doubleVector condensedTimes;
    intVector    condensedCycles;
};

STATE_API ostream& operator << (ostream &os, const DatabaseCorrelation &);

#endif
