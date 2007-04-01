#ifndef PLOT_H
#define PLOT_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: Plot
//
// Purpose:
//    This class is a plot element in a plot list.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Mar 30 12:25:01 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

class STATE_API Plot : public AttributeSubject
{
public:
    enum StateType
    {
        NewlyCreated,
        Pending,
        Completed,
        Error
    };

    Plot();
    Plot(const Plot &obj);
    virtual ~Plot();

    virtual void operator = (const Plot &obj);
    virtual bool operator == (const Plot &obj) const;
    virtual bool operator != (const Plot &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectPlotVar();
    void SelectDatabaseName();
    void SelectOperators();
    void SelectKeyframes();
    void SelectDatabaseKeyframes();

    // Property setting methods
    void SetStateType(StateType stateType_);
    void SetPlotType(int plotType_);
    void SetActiveFlag(bool activeFlag_);
    void SetHiddenFlag(bool hiddenFlag_);
    void SetExpandedFlag(bool expandedFlag_);
    void SetPlotVar(const std::string &plotVar_);
    void SetDatabaseName(const std::string &databaseName_);
    void SetOperators(const intVector &operators_);
    void SetActiveOperator(int activeOperator_);
    void SetId(int id_);
    void SetBeginFrame(int beginFrame_);
    void SetEndFrame(int endFrame_);
    void SetKeyframes(const intVector &keyframes_);
    void SetDatabaseKeyframes(const intVector &databaseKeyframes_);
    void SetIsFromSimulation(bool isFromSimulation_);

    // Property getting methods
    StateType         GetStateType() const;
    int               GetPlotType() const;
    bool              GetActiveFlag() const;
    bool              GetHiddenFlag() const;
    bool              GetExpandedFlag() const;
    const std::string &GetPlotVar() const;
          std::string &GetPlotVar();
    const std::string &GetDatabaseName() const;
          std::string &GetDatabaseName();
    const intVector   &GetOperators() const;
          intVector   &GetOperators();
    int               GetActiveOperator() const;
    int               GetId() const;
    int               GetBeginFrame() const;
    int               GetEndFrame() const;
    const intVector   &GetKeyframes() const;
          intVector   &GetKeyframes();
    const intVector   &GetDatabaseKeyframes() const;
          intVector   &GetDatabaseKeyframes();
    bool              GetIsFromSimulation() const;

    // Enum conversion functions
    static std::string StateType_ToString(StateType);
    static bool StateType_FromString(const std::string &, StateType &);
protected:
    static std::string StateType_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void AddOperator(int op);
    void ClearAllOperators();
    int  GetNumOperators() const;
    int  GetOperator(int i) const;
    void RemoveLastOperator();
private:
    int         stateType;
    int         plotType;
    bool        activeFlag;
    bool        hiddenFlag;
    bool        expandedFlag;
    std::string plotVar;
    std::string databaseName;
    intVector   operators;
    int         activeOperator;
    int         id;
    int         beginFrame;
    int         endFrame;
    intVector   keyframes;
    intVector   databaseKeyframes;
    bool        isFromSimulation;
};

#endif
