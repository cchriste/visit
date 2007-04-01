// ************************************************************************* //
//                         avtExpressionEvaluatorFilter.h                    //
// ************************************************************************* //

#ifndef AVT_EXPRESSION_EVALUATOR_FILTER_H
#define AVT_EXPRESSION_EVALUATOR_FILTER_H

#include <expression_exports.h>

#include <avtDatasetToDatasetFilter.h>
#include <avtQueryableSource.h>

#include <ExprPipelineState.h>


// ****************************************************************************
//  Class: avtExpressionEvaluatorFilter
//
//  Purpose:
//      This filter parses out expressions and turns them into executable
//      VTK networks.  It encapsulates the code that used to be in the
//      NetworkManager and the Viewer.
//
//  Programmer: Sean Ahern
//  Creation:   Thu Nov 21 15:15:07 PST 2002
//
// Modifications:
//   Brad Whitlock, Wed Aug 27 14:06:00 PST 2003
//   Made it use the right API.
//
//   Hank Childs, Mon Nov 17 16:47:33 PST 2003
//   Add ReleaseData.
//
//   Kathleen Bonnell, Thu Nov 13 08:39:40 PST 2003 
//   Added 'FindElementForPoint'.
//
//   Kathleen Bonnell, Mon Dec 22 14:39:30 PST 2003
//   Added GetDomainName.
//
//   Hank Childs, Thu Feb  5 17:11:06 PST 2004
//   Moved inlined constructor and destructor definitions to .C files
//   because certain compilers have problems with them.
//
//   Kathleen Bonnell, Tue May 25 16:16:25 PDT 2004 
//   Added QueryZoneCenter.
//
//   Kathleen Bonnell, Thu Jun 10 18:29:08 PDT 2004
//   Rename QueryZoneCenter to QueryCoords, added bool arg.
//
//   Kathleen Bonnell, Mon Jun 28 08:01:45 PDT 2004 
//   Added currentTimeState, ExamineSpecification. 
//
//   Kathleen Bonnell, Thu Dec 16 17:11:19 PST 2004 
//   Added another bool arg to QueryCoords. 
//
//   Hank Childs, Wed Dec 29 08:02:40 PST 2004
//   Added friend access to avtMacroExpressionFilter.  Also cache the
//   terminating source for updates.
//
// ****************************************************************************

class EXPRESSION_API avtExpressionEvaluatorFilter 
    : virtual public avtDatasetToDatasetFilter,
      virtual public avtQueryableSource
{
    friend class             avtMacroExpressionFilter;

  public:
                             avtExpressionEvaluatorFilter();
    virtual                 ~avtExpressionEvaluatorFilter();
    virtual const char*      GetType(void)
                                     { return "avtExpressionEvaluatorFilter";};

    virtual void             Query(PickAttributes *);
    virtual avtQueryableSource *
                             GetQueryableSource(void) { return this; };
    virtual void             ReleaseData(void);

    virtual bool             FindElementForPoint(const char *, const int, 
                                 const int, const char *, float[3], int &);
    virtual bool             QueryCoords(const std::string&, const int, 
                                 const int, const int, float[3], const bool,
                                 const bool);

    virtual void             GetDomainName(const std::string &, const int,
                                 const int , std::string &);

  protected:
    virtual void             PreExecute(void) {}
    virtual void             PostExecute(void) {}
    virtual void             Execute(void);
    virtual avtPipelineSpecification_p
                             PerformRestriction(avtPipelineSpecification_p);
    virtual int              AdditionalPipelineFilters(void);
    virtual void             ExamineSpecification(avtPipelineSpecification_p);

  protected:
    ExprPipelineState            pipelineState;
    avtPipelineSpecification_p   lastUsedSpec;
    avtSourceFromAVTDataset     *termsrc;

  private:
    int                          currentTimeState;
};


#endif


