// ************************************************************************* //
//                               avtMetaData.h                               //
// ************************************************************************* //

#ifndef AVT_META_DATA_H
#define AVT_META_DATA_H
#include <pipeline_exports.h>


#include <avtPipelineSpecification.h>

class     avtFacelist;
class     avtIntervalTree;
class     avtTerminatingSource;


// ****************************************************************************
//  Class: avtMetaData
//
//  Purpose:
//      This object is responsible for holding and possibly calculating meta-
//      data.  All of the meta-data is calculated in a lazy evaluation style -
//      only meta-data that is requested is calculated.  In addition, the
//      source is given an opportunity to give its copy of meta-data, in case
//      it is a source from a database and can go fetch some preprocessed 
//      meta-data.  The meta-data object currently only holds spatial and data
//      extents.
//
//  Programmer: Hank Childs
//  Creation:   August 7, 2000
//
//  Modifications:
//
//    Hank Childs, Tue Jun  5 08:05:29 PDT 2001
//    Changed interface to allow for a data specification to be considered
//    when getting materials or facelists.  Blew away outdated comments.
//
//    Jeremy Meredith, Thu Dec 13 16:00:29 PST 2001
//    Removed GetMaterial since it was never used.
//
// ****************************************************************************

class PIPELINE_API avtMetaData
{
  public:
                                 avtMetaData(avtTerminatingSource *);
    virtual                     ~avtMetaData() {;};

    avtIntervalTree             *GetDataExtents(void);
    avtIntervalTree             *GetSpatialExtents(void);

    avtFacelist                 *GetExternalFacelist(int);

  protected:
    avtTerminatingSource        *source;

    avtPipelineSpecification_p   GetPipelineSpecification(void);
    avtPipelineSpecification_p   GetPipelineSpecification(int);
};

#endif


