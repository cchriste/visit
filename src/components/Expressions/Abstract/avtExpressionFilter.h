// ************************************************************************* //
//                            avtExpressionFilter.h                          //
// ************************************************************************* //

#ifndef AVT_EXPRESSION_FILTER_H
#define AVT_EXPRESSION_FILTER_H

#include <expression_exports.h>

#include <string>

#include <avtStreamer.h>
#include <avtTypes.h>

class     vtkDataArray;
class     ArgsExpr;
class     ExprPipelineState;


// ****************************************************************************
//  Class: avtExpressionFilter
//
//  Purpose:
//      This is a base class that lets derived types not worry about how to set
//      up a derived variable.
//
//  Programmer: Hank Childs
//  Creation:   June 7, 2002
//
//  Modifications:
//
//    Sean Ahern, Fri Jun 13 11:18:07 PDT 2003
//    Added the virtual function NumVariableArguments that lets an
//    expression declare how many of its arguments are variables.
//
//    Hank Childs, Tue Aug 12 10:34:02 PDT 2003
//    Store the domain label and index when deriving a variable.  Certain
//    filters need this information.
//  
//    Hank Childs, Fri Oct 24 14:46:20 PDT 2003
//    Made Evaluator be a friend class so that it could call perform
//    restriction on the expressions it contains.
//
//    Hank Childs, Wed Dec 10 09:40:18 PST 2003
//    Add support for recentering a variable.  It is at this level so all
//    derived types can use the same routine.
// 
//    Kathleen Bonnell, Mon Jun 28 07:48:55 PDT 2004 
//    Add currentTimeState, ExamineSpecification.
//
//    Hank Childs, Mon Dec 27 10:13:44 PST 2004
//    Separate out parts that are related to streaming.
//
//    Hank Childs, Fri Aug  5 16:40:12 PDT 2005
//    Added GetVariableType.
//
//    Hank Childs, Mon Aug 29 14:44:20 PDT 2005
//    Added SetExpressionAtts.
//
// ****************************************************************************

class EXPRESSION_API avtExpressionFilter : virtual public 
                                                     avtDatasetToDatasetFilter
{
    friend class             avtExpressionEvaluatorFilter;

  public:
                             avtExpressionFilter();
    virtual                 ~avtExpressionFilter();

    void                     SetOutputVariableName(const char *);
    virtual void             AddInputVariableName(const char *var)
                                {SetActiveVariable(var);}

    virtual void             ProcessArguments(ArgsExpr *, ExprPipelineState *);
    virtual int              NumVariableArguments() = 0;

  protected:
    char                    *outputVariableName;
    int                      currentTimeState;

    virtual bool             IsPointVariable();

    virtual void             PreExecute(void);
    virtual void             PostExecute(void);
    virtual void             RefashionDataObjectInfo(void);
    void                     SetExpressionAttributes(const avtDataAttributes &,
                                                     avtDataAttributes &);
    virtual avtPipelineSpecification_p
                             PerformRestriction(avtPipelineSpecification_p);
    virtual void             ExamineSpecification(avtPipelineSpecification_p);

    virtual int              GetVariableDimension() { return 1; };
    virtual avtVarType       GetVariableType() { return AVT_UNKNOWN_TYPE; };

    vtkDataArray            *Recenter(vtkDataSet*, vtkDataArray*,avtCentering);
    void                     UpdateExtents(avtDataTree_p);
};


#endif


