// ************************************************************************* //
//                           avtFeatureEdgesFilter.h                         //
// ************************************************************************* //

#ifndef AVT_FEATURE_EDGES_FILTER_H
#define AVT_FEATURE_EDGES_FILTER_H

#include <pipeline_exports.h>

#include <avtStreamer.h>


// ****************************************************************************
//  Class: avtFeatureEdgesFilter
//
//  Purpose:
//      Extracts feature edges of a dataset.
//
//  Programmer: Jeremy Meredith
//  Creation:   March 12, 2002
//
//  Modifications:
//
//    Hank Childs, Tue Sep 10 15:31:39 PDT 2002
//    Improve the way memory is handled.  Also removed stuff for making this
//    be a plugin filter and inherited from avtStreamer.
//
// ****************************************************************************

class PIPELINE_API avtFeatureEdgesFilter : public avtStreamer
{
  public:
                             avtFeatureEdgesFilter() {;};
    virtual                 ~avtFeatureEdgesFilter() {;};

    virtual const char      *GetType(void) { return "avtFeatureEdgesFilter"; };
    virtual const char      *GetDescription(void)
                                 { return "Extracting feature edges"; };

  protected:
    virtual vtkDataSet      *ExecuteData(vtkDataSet *, int, std::string);
    virtual void             RefashionDataObjectInfo(void);
};


#endif


