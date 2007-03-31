// ************************************************************************* //
//                            avtSummationQuery.h                            //
// ************************************************************************* //

#ifndef AVT_SUMMATION_QUERY_H
#define AVT_SUMMATION_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>

#include <string>

class vtkDataSet;

// ****************************************************************************
//  Class: avtSummationQuery
//
//  Purpose:
//      This query sums all of the values for a variable.
//
//  Notes:
//    Taken mostly from Hank Childs' avtSummationFilter and reworked to
//    fit into the Query hierarchy.  
//
//  Programmer: Kathleen Bonnell 
//  Creation:   September 30, 2002 
//
//  Modifications:
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002
//    Add domain to Execute arguments.
//
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002
//    Add unitsAppend. 
//
// ****************************************************************************

class QUERY_API avtSummationQuery : public avtDatasetQuery
{
  public:
                                    avtSummationQuery();
    virtual                        ~avtSummationQuery() {;};

    virtual void                    SetVariableName(std::string &);
    void                            SetSumType(std::string &);
    void                            SetUnitsAppend(std::string &);

    virtual const char             *GetType(void)
                                             { return "avtSummationQuery"; };
    virtual const char             *GetDescription(void)
                                             { return descriptionBuffer; };

    void                            SumGhostValues(bool);
    void                            SumOnlyPositiveValues(bool);

  protected:
    double                          sum;
    std::string                     variableName;
    std::string                     sumType;
    std::string                     unitsAppend;
    bool                            sumGhostValues;
    bool                            sumOnlyPositiveValues;
    char                            descriptionBuffer[1024];

    virtual void                    Execute(vtkDataSet *, const int);
    virtual void                    PreExecute(void);
    virtual void                    PostExecute(void);
};


#endif


