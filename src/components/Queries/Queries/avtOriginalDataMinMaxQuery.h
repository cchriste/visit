// ************************************************************************* //
//                     avtOriginalDataMinMaxQuery.h                          //
// ************************************************************************* //

#ifndef AVT_ORIGINAL_MINMAX_QUERY_H
#define AVT_ORIGINAL_MINMAX_QUERY_H
#include <query_exports.h>

#include <avtMinMaxQuery.h>

class avtExpressionEvaluatorFilter;


// ****************************************************************************
//  Class: avtOriginalDataMinMaxQuery
//
//  Purpose:
//    A query that retrieves the min/max of a var from the original data. 
//
//  Programmer: Kathleen Bonnell
//  Creation:   February 10, 2004 
//
//  Modifications:
//
// ****************************************************************************

class QUERY_API avtOriginalDataMinMaxQuery : public avtMinMaxQuery
{
  public:
                                  avtOriginalDataMinMaxQuery();
    virtual                      ~avtOriginalDataMinMaxQuery();

    virtual bool                  OriginalData(void) { return true; };

  protected:
    virtual avtDataObject_p       ApplyFilters(avtDataObject_p);   
    virtual int                   GetNFilters() { return 1; };

  private:
    avtExpressionEvaluatorFilter *eef;
};

#endif
