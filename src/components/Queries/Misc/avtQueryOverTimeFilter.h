// ************************************************************************* //
//                       avtQueryOverTimeFilter.h                            //
// ************************************************************************* //

#ifndef AVT_QUERYOVERTIME_FILTER_H
#define AVT_QUERYOVERTIME_FILTER_H

#include <query_exports.h>

#include <avtDatasetToDatasetFilter.h>
#include <avtTimeLoopFilter.h>

#include <QueryOverTimeAttributes.h>
#include <SILRestrictionAttributes.h>

class vtkPolyData;

// ****************************************************************************
//  Class: avtQueryOverTimeFilter
//
//  Purpose:
//    Performs a query over time. 
//
//  Programmer: Kathleen Bonnell 
//  Creation:   March 19, 2004
//
//  Modifications:
//    Brad Whitlock, Wed Apr 14 14:56:45 PST 2004
//    Fixed for Windows.
//
//    Kathleen Bonnell, Tue May  4 14:21:37 PDT 2004
//    Removed SilUseSet in favor of SILRestrictionAttributes. 
//
//    Kathleen Bonnell, Thu Jan  6 11:12:35 PST 2005 
//    Added inheritance from avtTimeLoopFilter, which handles the stepping
//    through time.  Removed PostExecute method.  Added CreateFinalOutput,
//    ExecutionSuccessful (required by new inheritance).  Added qRes, times
//    (used to be stored by avtDataSetQuery::PerformQueryInTime).  Added
//    sucess and finalOutputCreated.
//
//    Kathleen Bonnell, Tue Nov  8 10:45:43 PST 2005
//    Added CreatePolys method, and members useTimeForXAxis, nResultsToStore. 
//
// ****************************************************************************

class QUERY_API avtQueryOverTimeFilter : public avtTimeLoopFilter,
                                         public avtDatasetToDatasetFilter
{
  public:
                          avtQueryOverTimeFilter(const AttributeGroup*);
    virtual              ~avtQueryOverTimeFilter();

    static avtFilter     *Create(const AttributeGroup*);

    virtual const char   *GetType(void)  { return "avtQueryOverTimeFilter"; };
    virtual const char   *GetDescription(void) { return "Querying over Time"; };

    void                  SetSILAtts(const SILRestrictionAttributes *silAtts);

  protected:
    QueryOverTimeAttributes   atts;
    SILRestrictionAttributes  querySILAtts;
    doubleVector          qRes;
    doubleVector          times;
    bool                  success;
    bool                  finalOutputCreated;

    bool                  useTimeForXAxis;
    int                   nResultsToStore;

    virtual void          Execute(void);
    virtual void          RefashionDataObjectInfo(void);

    virtual int           AdditionalPipelineFilters(void) { return 1; };

    virtual void          CreateFinalOutput(void);
    virtual bool          ExecutionSuccessful(void) { return success; } ;
    vtkPolyData          *CreatePolys(const doubleVector &, 
                                      const doubleVector &);
};


#endif


