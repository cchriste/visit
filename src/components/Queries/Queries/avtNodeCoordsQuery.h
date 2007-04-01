// ************************************************************************* //
//                          avtNodeCoordsQuery.h                             //
// ************************************************************************* //

#ifndef AVT_NODECOORDS_QUERY_H
#define AVT_NODECOORDS_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>


// ****************************************************************************
//  Class: avtNodeCoordsQuery
//
//  Purpose:
//      This is a dataset query that returns the coordintes of a node. 
//
//  Programmer: Kathleen Bonnell 
//  Creation:   June 10, 2004 
//
//  Modifications:
//    Kathleen Bonnell, Thu Dec 16 17:16:33 PST 2004
//    Add FindLocalCoord and FindGlobalCoord.
//
// ****************************************************************************

class QUERY_API avtNodeCoordsQuery : public avtDatasetQuery
{
  public:
                             avtNodeCoordsQuery();
    virtual                 ~avtNodeCoordsQuery(); 

    virtual const char      *GetType(void) { return "avtNodeCoordsQuery"; };
    virtual const char      *GetDescription(void) 
                                  { return "Getting node coords."; };

    virtual void             PerformQuery(QueryAttributes *);
    virtual bool             OriginalData(void) { return true; };

  protected:

    virtual void             Execute(vtkDataSet*, const int){;};
    bool                     FindLocalCoord(float[3]);
    bool                     FindGlobalCoord(float[3]);
};

#endif

